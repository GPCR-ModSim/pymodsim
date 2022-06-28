import os

# This is the folder where pymemdyn git repo has been deployed,
# or to be more specific, where the settings.py file is located at.
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The following variable can be given as a full path to the place where
# the templates reside or it can also use a relative path to the previously
# defined ROOT_DIR.
TEMPLATES_DIR = os.path.join(ROOT_DIR, "templates")


# AlphaFold settings
# See https://github.com/deepmind/alphafold/ for more information on each of
# the settings.

# Define the path of alphafold2.sif:
#ALPHAFOLD_SIF = {your_ALPHAFOLD_SIF}

# Define the path of the genetic databases:
#ALPHAFOLD_DATA_PATH = {your_ALPHAFOLD_DATA_PATH}      

# Define the path of the AlphaFold parameters:
#ALPHAFOLD_MODELS = {your_ALPHAFOLD_MODELS}

# Define the path of the genetic databases in relation to the genetic database
# path (defined as /data):
#BFD_DATABASE_PATH = {your_BFD_DATABASE_PATH}
#MGNIFY_DATABASE_PATH = {your_MGNIFY_DATABASE_PATH}
#PDB70_DATABASE_PATH = {your_PDB70_DATABASE_PATH}
#TEMPLATE_MMCIF_DIR = {your_TEMPLATE_MMCIF_DIR}
#OBSOLETE_PDBS_PATH = {your_OBSOLETE_PDBS_PATH}
#UNICLUST30_DATABASE_PATH = {your_UNICLUST30_DATABASE_PATH}
#UNIREF90_DATABASE_PATH = {your_UNIREF90_DATABASE_PATH}

# ALPHAFOLD_PATHs examples:
ALPHAFOLD_SIF = "/home/apps/apps/alphafold/alphafold2.sif"
ALPHAFOLD_DATA_PATH = "/home/apps/apps/alphafold/data_2021-12-16"   
ALPHAFOLD_MODELS = "/home/apps/apps/alphafold/data_2021-12-16/params"
BFD_DATABASE_PATH = "/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
MGNIFY_DATABASE_PATH = "/data/mgnify/mgy_clusters_2018_12.fa"
PDB70_DATABASE_PATH = "/data/pdb70/pdb70"
TEMPLATE_MMCIF_DIR = "/data/pdb_mmcif/mmcif_files"
OBSOLETE_PDBS_PATH = "/data/pdb_mmcif/obsolete.dat"
UNICLUST30_DATABASE_PATH = "/data/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
UNIREF90_DATABASE_PATH = "/data/uniref90/uniref90.fasta"


# PPM settings

# Define the path of PPM 3.0
#PPM_PATH = {your_PPM_PATH}

# PPM_PATH example:
PPM_PATH = "/home/rbroek/apps/PPM_3.0"
