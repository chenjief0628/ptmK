# ptmK
A pipeline is used to add side chain modifications to lysine

## Overview
Post-translational modifications (PTMs) expand protein functionality and diversity, which leads to increased proteome complexity and lysine side-chain modifications being one of most prevalent amino acids. Molecular docking and molecular dynamics simulations are widely used to identify the key amino acids from among numerous modified amino acids. However, existing tools lack accessible solutions for high-throughput generating protein structures with different post-translational side-chain modifications. To address this gap, we have developed ptmK, a user-friendly pipeline specifically designed to generate structures of lysine with post-translational side-chain modifications.

# Installation
We store the public release versions of ptmK on GitHub, a site that provides code development with version control and issue tracking through the use of git. We will not describe the use of git in general, as you will not need more than very basic features. Below we outline the few commands needed on a windows system; please refer to general git descriptions and tutorials to suit your system. To get the code, you clone or download the repository. We recommend cloning, as it allows you to easily update the code when new versions are released.
## Installation of ptmK
bioconda is required
```
git clone https://github.com/chenjief0628/ptmk/
conda install python=3.12.3
pip install pandas
pip install biopython
conda install -c conda-forge pdbfixer openmm -y
conda install freesasa
```
# Useful Example
The GP6D structure 2bh9.pdb was used for test. 
ptmK consists of three main stages: 
1. Preprocessing the structure by removing non-amino acid atoms and completing any missing amino acid residues.
2. Calculating the modification potential for each lysine residue, which includes assessing the solvent-accessible surface area of each lysine and its distance to the active site; if the location of the active site is unknown, lysine residues can be selected randomly. The final result saved as cal_distance_sasa.csv in the current directory.
3. Applying the specified modifications to the designated amino acids, with the modification structures stored in the motif folder and the final modified results saved as motif.result.pdb in the current directory.

## Notice
Some PDB files do not contain full-length amino acid sequences. After the cleaning process, the residue numbering will be reset. Therefore, it is crucial to carefully verify the sequence position of the modification site relative to the actual sequence position in the PDB file. If the lysine residue that requires modification has been identified, the third step can be initiated directly.
```
python clean.py --pdb_path ./G6PD/2bh9.pdb --out_path ./G6PD/2bh9_A_chain.pdb
python cla_sasa.py ./G6PD/2bh9_A_chain.pdb --residue_ids 170,149,243
python motif.py --pdb_path ./G6PD/2bh9_A_chain.pdb --motif_path ./motif/KAC.pdb --nums 19
```
