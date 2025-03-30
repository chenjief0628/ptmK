import os
import re
import uuid
import shutil
import time
import json
import sqlite3
import threading
from fastapi import FastAPI, File, UploadFile, Form, BackgroundTasks, Request, HTTPException
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware
from fastapi.middleware.cors import CORSMiddleware
#front
from fastapi.staticfiles import StaticFiles
from fastapi.responses import RedirectResponse

from Bio import PDB
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import numpy as np
import pandas as pd
import freesasa

app = FastAPI()
origins = [
    "http://127.0.0.1:8000",
    "http://localhost:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,          
    allow_credentials=True,         
    allow_methods=["*"],            
    allow_headers=["*"],            
)

# savedir
UPLOAD_DIR = "uploads"
PROJECT_DIR = "projects"
RESULT_DIR = "results"
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(PROJECT_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)

######################################
# 1 RateLimit
######################################
class RateLimitMiddleware(BaseHTTPMiddleware):
    def __init__(self, app, window_seconds: int = 60, max_total_upload_size: int = 10 * 1024 * 1024, max_tasks: int = 5):
        super().__init__(app)
        self.window_seconds = window_seconds
        self.max_total_upload_size = max_total_upload_size
        self.max_tasks = max_tasks
        self.ip_tasks = {}

    async def dispatch(self, request: Request, call_next):
        ip = request.client.host
        current_time = time.time()
        
        if request.url.path == "/upload" and request.method.upper() == "POST":
            try:
                upload_size = int(request.headers.get("content-length", 0))
            except ValueError:
                upload_size = 0

            task_records = self.ip_tasks.get(ip, [])
            task_records = [record for record in task_records if current_time - record[0] < self.window_seconds]
            
            total_tasks = len(task_records)
            total_uploaded = sum(record[1] for record in task_records)
            
            if total_tasks >= self.max_tasks:
                return JSONResponse(status_code=429, content={"detail": "任务提交过多，请稍后重试"})
            if total_uploaded + upload_size > self.max_total_upload_size:
                return JSONResponse(status_code=429, content={"detail": "上传文件总量超出限制，请稍后重试"})
            
            task_records.append((current_time, upload_size))
            self.ip_tasks[ip] = task_records

        response = await call_next(request)
        return response

app.add_middleware(RateLimitMiddleware, window_seconds=60, max_total_upload_size=10 * 1024 * 1024, max_tasks=5)

######################################
# 2 secure filename
######################################
def secure_filename(filename: str) -> str:
    filename = os.path.basename(filename)
    filename = re.sub(r"[^A-Za-z0-9_.-]", "", filename)
    return filename

######################################
# 3 SQLite
######################################
DB_PATH = "proj.db"
db_lock = threading.Lock()

def init_db():
    need_recreate = False
    try:
        conn = sqlite3.connect(DB_PATH, check_same_thread=False, timeout=10)
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA foreign_keys=ON;")
        cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='projects';")
        if not cursor.fetchone():
            need_recreate = True
    except Exception as e:
        need_recreate = True

    if need_recreate:
        try:
            if os.path.exists(DB_PATH):
                conn.close()
                os.remove(DB_PATH)
        except Exception:
            pass
        conn = sqlite3.connect(DB_PATH, check_same_thread=False, timeout=10)
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA foreign_keys=ON;")
        conn.execute('''CREATE TABLE projects (
                        project_id TEXT PRIMARY KEY,
                        status TEXT,
                        progress TEXT,
                        error TEXT,
                        files TEXT)''')
        conn.commit()
    return conn

db_conn = init_db()

def insert_project(project_id, status, progress, error=None, files=None):
    with db_lock:
        db_conn.execute(
            "INSERT INTO projects (project_id, status, progress, error, files) VALUES (?, ?, ?, ?, ?)",
            (project_id, status, progress, error, json.dumps(files) if files else None)
        )
        db_conn.commit()

def update_project(project_id, status=None, progress=None, error=None, files=None):
    with db_lock:
        fields = []
        params = []
        if status is not None:
            fields.append("status = ?")
            params.append(status)
        if progress is not None:
            fields.append("progress = ?")
            params.append(progress)
        if error is not None:
            fields.append("error = ?")
            params.append(error)
        if files is not None:
            fields.append("files = ?")
            params.append(json.dumps(files))
        if not fields:
            return
        params.append(project_id)
        db_conn.execute("UPDATE projects SET " + ", ".join(fields) + " WHERE project_id = ?", params)
        db_conn.commit()

def get_project(project_id):
    with db_lock:
        cursor = db_conn.execute("SELECT project_id, status, progress, error, files FROM projects WHERE project_id = ?", (project_id,))
        row = cursor.fetchone()
        if row:
            return {
                "project_id": row[0],
                "status": row[1],
                "progress": row[2],
                "error": row[3],
                "files": json.loads(row[4]) if row[4] else {}
            }
        else:
            return None

######################################
# 4 Process
######################################

#clean.py modified
import os
from Bio import PDB
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import argparse

class NonAminoAcidRemover(PDB.Select):
    def accept_residue(self, residue):
        # 检查残基是否是氨基酸
        if residue.get_id()[0] == ' ':
            return 1  # 是氨基酸，保留
        return 0  # 不是氨基酸，删除

def process_ent_file(ent_file_path, output_pdb_path):
    # 1. 使用 BioPython 读取 .ent 文件，提取 A 链
    parser = PDB.PPBuilder()
    
    # 读取 .ent 文件并解析
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', ent_file_path)
    model = structure[0]
    chains = []
    for model in structure:
        for chain in model:
            id = chain.get_id()
            chains.append(id)
    if len(chains) > 1:
        print("chain id are " + ",".join(chains) + " and first chain was selected.")
    else:
        print("chan id is " + chains[0])
    a_chain = model[chains[0]]
    
    with open("temp.pdb",'w') as f:
        io = PDB.PDBIO()
        io.set_structure(a_chain)
        io.save(f,NonAminoAcidRemover())
        
    fixer = PDBFixer(filename= "temp.pdb")
    
    # 修复 A 链的缺失原子
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # 添加氢原子，使用7.0的pH值
    
    # 3. 将修复后的结构保存为 PDB 文件
    # 保存修复后的 PDB 文件
    with open(output_pdb_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"A 链已修复并保存为 {output_pdb_path}")

def run_clean(input_file: str, output_file: str):
    process_ent_file(input_file, output_file)


#cla_sasa.py modified
from Bio.PDB import PDBParser, is_aa
import numpy as np
import os
import pandas as pd
from Bio.PDB import DSSP
import re
import freesasa
import argparse
def get_lys_residues(model):
    lys_residues = {}
    for chain in model:
        for residue in chain:
            if residue.get_resname() == 'LYS':
                res_id = residue.get_id()
                key = (chain.id, res_id[1], res_id[2].strip())
                if 'CA' in residue:
                    lys_residues[key] = residue['CA'].get_coord()
                else:
                    print("There is no CA atom in " + chain.id + " " + res_id[1])
    return lys_residues


def calculate_centroid(model, residue_ids):
    ca_coords = []
    for chain in model:
        chain_id = chain.get_id()
    
    for residue_id in residue_ids:
        residue = model[chain_id][(' ', int(residue_id) -1 , ' ')]
        # 获取 CA 原子
        ca = residue['CA']
        coord = ca.get_coord()
        ca_coords.append(coord)
        #print(ca_coords)
    
    centroid = np.mean(ca_coords, axis=0)
    return centroid

def cal_dist(dict_k, zx):
    list_d = []
    for key in dict_k:
        c1 = dict_k[key]
        c2 = zx
        d = np.linalg.norm(c1 - c2)
        list_d.append(d)
    return list_d

def cal_sasa(pdb_file,model):
    aim_list = []
    for chain in model:
        for residue in chain:
            if residue.get_resname() == 'LYS':
                res_id = residue.get_id()
                aim = res_id[1]
                aim_list.append(aim)
            
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    result = result.residueAreas()['A']
    vlist = []
    for aim in aim_list:
        v = result[str(aim)].total
        vlist.append(v)
    return vlist

def run_cla_sasa(input_file: str, output_csv: str, residue_ids: str):
    
    pdb_file = input_file
    residue_ids = residue_ids.split(',')
    #residue_ids_list = residue_ids.split(',')
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]
    ds = calculate_centroid(model, residue_ids)
    dict_k = get_lys_residues(model)
    d = cal_dist(dict_k,ds)
    vl = cal_sasa(pdb_file, model)
    list_k = []
    for k in dict_k:
        num = k[1]
        list_k.append(num)
    score = [1/a * 0.7 + b/max(vl) * 0.3 for a, b in zip(d, vl)]
    #print(score)
    data = {'K num': list_k, 'Activate distance': d, 'SASA area':vl, 'score' : score}
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"计算完成，结果已保存到 {output_csv}")

#motif.py modified
import os
from Bio import PDB
import numpy as np
import pandas as pd
import argparse
from Bio.PDB import PDBParser, PDBIO

def get_residue(pdb_path):
    """
    Get the residue coordinates list for pdb file
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_path)
    model = structure[0]
    for chain in model:
        for residue in chain:
            atom_coordinates = {}
            for atom in residue:
                atom1 = atom.get_name()
                atom_coordinates[atom1] = atom.get_coord()
    return atom_coordinates

def center_coordinates(coordinates):
    """
    Get the center coordinates of the residue
    """
    center = np.mean(coordinates, axis=0)
    return coordinates - center


def move_coord(motif_coords, pdb_path, nums):
    """
    Align the motif coordinates to the pdb coordinates
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_path)
    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.get_id()[1] == nums:
                atom_coordinates = {}
                for atom in residue:
                    atom1 = atom.get_name()
                    atom_coordinates[atom1] = atom.get_coord()
    pdb_coords = atom_coordinates
    overlap = ['N', 'CA', 'C', 'O']
    coords_dict1 = [pdb_coords[atom] for atom in overlap if atom in pdb_coords]
    coords_dict2 = [motif_coords[atom] for atom in overlap if atom in motif_coords]
    coords_dict1 = np.array(coords_dict1)
    coords_dict2 = np.array(coords_dict2)
    centroid1 = np.mean(coords_dict1, axis=0)
    centroid2 = np.mean(coords_dict2, axis=0)
    coords_dict1_translated = coords_dict1 - centroid1
    coords_dict2_translated = coords_dict2 - centroid2
    H = np.dot(coords_dict2_translated.T, coords_dict1_translated)
    U, _, Vt = np.linalg.svd(H)
    rotation_matrix = np.dot(Vt.T, U.T)
    transformed_coords = {}
    for atom, coord in motif_coords.items():
        # 平移
        coord_translated = coord - centroid2
        # 旋转
        coord_rotated = np.dot(rotation_matrix, coord_translated.T) + centroid1
        transformed_coords[atom] = coord_rotated

    return transformed_coords

def save_to_pdb(coords_dict, pdb_path, motif_path, nums,output_file):
    """
    Save the transformed coordinates to a new pdb file
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', motif_path)
    model = structure[0]
    for chain in model:
        for residue in chain:
            residue.id = (' ', nums, ' ')
            for atom_name, new_coord in coords_dict.items():
                atom = residue[atom_name]
                atom.coord = np.array(new_coord)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save('temp.pdb')
    structure1 = parser.get_structure('PDB', pdb_path)
    model1 = structure1[0]
    for chain in model1:
        residues = list(chain)
        residue_to_delete = residues[nums - 1]
        chain.detach_child(residue_to_delete.get_id())
        residues = list(chain)
        residues.insert(nums -1, residue)
    
    aaa = list(chain)
    for i in aaa:
        chain.detach_child(i.get_id())

    for i in residues:
        #print(i)
        chain.add(i)
    io = PDB.PDBIO()
    io.set_structure(structure1)
    io.save(output_file)
    print(f"motif 对齐完成，生成结果文件 {output_file}")

def run_motif(input_file: str, motif_file: str, output_file: str, residue_num: int):
    motif_path=motif_file
    pdb_path=input_file
    nums=residue_num
    motif_coords = get_residue(motif_path)
    transformed = move_coord(motif_coords, pdb_path, nums)
    #print(transformed)
    save_to_pdb(transformed, pdb_path, motif_path,nums,output_file)


# run process
def process_project(project_id: str, uploaded_filepath: str, residue_ids: str, motif_file: str, residue_num: int):
    proj_dir = os.path.join(PROJECT_DIR, project_id)
    os.makedirs(proj_dir, exist_ok=True)
    try:
        # clean
        clean_output = os.path.join(proj_dir, "cleaned.pdb")
        run_clean(uploaded_filepath, clean_output)
        update_project(project_id, progress="cleaned")
        
        # calc_sasa
        csv_output = os.path.join(proj_dir, "cal_distance_sasa.csv")
        run_cla_sasa(clean_output, csv_output, residue_ids)
        update_project(project_id, progress="calculated")
        
        # motif
        motif_output = os.path.join(proj_dir, "motif.result.pdb")
        run_motif(clean_output, motif_file, motif_output, residue_num)
        update_project(project_id, progress="modified")
        
        # update
        update_project(project_id, status="completed", files={
            "cleaned_pdb": clean_output,
            "cal_distance_sasa_csv": csv_output,
            "motif_result_pdb": motif_output
        })
    except Exception as e:
        update_project(project_id, status="failed", error=str(e))

######################################
# 5 API
######################################
@app.post("/upload")
async def upload_file(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    residue_ids: str = Form(...),
    motif_file: UploadFile = File(...),
    residue_num: int = Form(...)
):
    filename = secure_filename(file.filename)
    motif_filename = secure_filename(motif_file.filename)
    
    project_id = str(uuid.uuid4())
    upload_path = os.path.join(UPLOAD_DIR, project_id)
    os.makedirs(upload_path, exist_ok=True)
    file_location = os.path.join(upload_path, filename)
    motif_location = os.path.join(upload_path, motif_filename)
    
    with open(file_location, "wb") as f:
        f.write(await file.read())
    with open(motif_location, "wb") as f:
        f.write(await motif_file.read())
    
    insert_project(project_id, "processing", "init", files={})
    
    background_tasks.add_task(process_project, project_id, file_location, residue_ids, motif_location, residue_num)
    
    return {"project_id": project_id, "message": "Proj Start"}

@app.get("/project/{project_id}")
def get_project_status(project_id: str):
    project = get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Proj Not Found")
    return project

@app.get("/download/{project_id}/{file_key}")
def download_file(project_id: str, file_key: str):
    project = get_project(project_id)
    if not project or project["status"] != "completed":
        raise HTTPException(status_code=404, detail="Proj Not Found or Unfinished")
    files = project["files"]
    if file_key not in files:
        raise HTTPException(status_code=404, detail="File Not Found")
    return FileResponse(files[file_key], filename=os.path.basename(files[file_key]))

#front
app.mount("/front", StaticFiles(directory="/app/front", html=True), name="front")
@app.get("/", response_class=HTMLResponse)
def read_index():
    return RedirectResponse(url="/front/index.html")

######################################
# 6 server
######################################
@app.on_event("shutdown")
def shutdown_event():
    with db_lock:
        db_conn.commit()
        db_conn.close()
    print("DB closed, server shutting down")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
