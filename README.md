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
pip install -r requirements.txt
```
if python packages pdbfixer openmm and freesasa can not be installed , you can install it using conda
```
conda install -c conda-forge pdbfixer openmm -y
conda install freesasa
```
# Useful Example

```
python ./code/clean.py --pdb_path ./G6PD/2bh9.pdb --out_path ./G6PD/2bh9_A_chain.pdb

```
