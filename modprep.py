import modeller
from modeller import automodel
import sys


def run_modeller(**kwargs):
    """
    Refine protein structure using MODELLER
    """
    default_stdout = sys.stdout
    
    sys.stdout = open("pymodsim.log", "a")
    
    if kwargs["mode"] == "run":
        modeller.log.minimal()
        env = modeller.Environ()
        env.io.atom_files_directory = ['.', '../atom_files']
        a = automodel.AutoModel(env,
                               alnfile  = kwargs["alnfile"],  
                               knowns   = kwargs["knowns"],     
                               sequence = 'refined')
        a.starting_model= 1              
        a.ending_model  = 1              
        a.md_level = None                 
        
        a.make()
        
    sys.stdout = default_stdout
    