##########################################################################
#                 Homology                                               #
##########################################################################

class Homology(object):
    def __init__(self, **kwargs):
        self.steps = ["clean_fasta", "set_stage_init", "alphafold", "set_stage_init2", "set_end"]
 
        self.recipe = \
            {"clean_fasta":     {"command":    "clean_fasta",  # 1
                                 "options":   {"seq":  ""}},
                
             "set_stage_init":  {"command":    "set_stage_init",  # 2
                                 "options":   {"src_dir":    "",
                                               "src_files":  "",
                                               "tgt_dir":    "sequences"}},

             "alphafold":       {"alphafold":  "alphafold",  # 3
                                 "options":   {"seq": "",
                                               "max_template_date": "2021-11-01"}},
             
             "set_stage_init2": {"command":    "set_stage_init",  # 4
                                 "options":   {"src_dir":    "",
                                               "src_files": ["ranked_0.pdb"],
                                               "tgt_dir":    "."}},
             
             "set_end":         {"command":    "set_stage_init",  # 5
                                 "options":   {"src_dir":    "",
                                               "src_files": ["ranked_0.pdb"],
                                               "tgt_dir":    "finalOutput"}}}
                                 
        self.breaks = \
            {"clean_fasta":     {"seq":       "sequence"},
             "set_stage_init":  {"src_files": "sequence"},
             "alphafold":       {"seq":       "sequence"},
             "set_stage_init2": {"src_dir":   "sequence_base"}}          


##########################################################################
#                 ModPrep                                                #
##########################################################################

class ModPrep(object):
    def __init__(self, **kwargs):
        self.steps = ["clean_fasta", "plot_conf", "make_pir", "run_modeller", "set_end"]
        
        self.recipe = \
           {"clean_fasta":     {"command":  "clean_fasta",  # 1
                                "options": {"seq":  ""}},
            
            "plot_conf":       {"command":  "plot_conf",  # 2
                                "options": {"pdb": "",
                                             "tgt": "plot_confidence.txt"}},
                
            "make_pir":        {"command":  "make_pir",  # 3
                                "options": {"seq": "",
                                            "pdb": "",
                                            "tgt1": "alignment.pir",
                                            "tgt2": "refinement.txt",
                                            "tgt3": "alignment.txt"}},
                         
            "run_modeller":    {"command":   "run_modeller",  # 4
                                "options":  {"alnfile": "alignment.pir",
                                             "knowns":  "",
                                             "mode":    ""}},
                        
            "set_end":         {"command":    "set_stage_init",  # 5
                                "options":   {"src_dir":    "",
                                              "src_files": ["plot_confidence.txt",
                                                            "alignment.txt",
                                                            "refinement.txt",
                                                            "refined.B99990001.pdb"],
                                              "tgt_dir":    "finalOutput"}}}                                
            
        self.breaks = \
            {"clean_fasta":      {"seq": "sequence"},
             "plot_conf":       {"pdb": "pdb_2"},
             "make_pir":        {"seq": "sequence",
                                 "pdb": "pdb_2"},
             "run_modeller":    {"knowns": "pdb_2",
                                 "mode": "mode"}}


##########################################################################
#                 Orientation                                            #
##########################################################################

class Orientation(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "get_protein", "make_inp", "ppm", "clean_pdb", 
                      "superimpose", "set_end1", "set_end2"]
        
        self.recipe = \
            {"set_stage_init":  {"command":  "set_stage_init",  # 1
                                 "options": {"repo_files": ["res.lib"],
                                             "tgt_dir":     "."}},
             
             "get_protein":     {"command":  "get_protein",  # 2
                                 "options": {"pdb": "",
                                             "pdb_backup": "",
                                             "tgt": "protein_stripped.pdb"}},

             "make_inp":        {"command":  "make_inp",  # 3
                                 "options": {"pdb": "protein_stripped.pdb",
                                             "tgt": "2membranes.inp"}},
           
             "ppm":             {"ppm":      "ppm",  # 4
                                 "options": {"inp": "2membranes.inp"}},
             
             "clean_pdb":       {"command":  "clean_pdb",  # 5
                                 "options": {"src": "protein_strippedout.pdb",
                                             "tgt": "protein_aligned.pdb"}},
             
             "superimpose":     {"command":  "superimpose",  # 6
                                 "options": {"pdb": "",
                                             "pdb_backup": "",
                                             "src": "protein_aligned.pdb",
                                             "tgt": "homology.pdb"}},

             "set_end1":         {"command":   "set_stage_init",  # 7
                                 "options":   {"src_dir":    "",
                                               "src_files":  "",
                                               "tgt_dir":    "finalOutput"}},
             
             "set_end2":         {"command":    "set_stage_init",  # 8
                                 "options":   {"src_dir":    "",
                                               "src_files":  ["protein_stripped.pdb",
                                                              "protein_aligned.pdb",
                                                              "homology.pdb"],
                                               "tgt_dir":    "finalOutput"}}}
            
        self.breaks = \
            {"get_protein":     {"pdb": "pdb_3",
                                 "pdb_backup": "pdb_2"},
             "superimpose":     {"pdb": "pdb_3",
                                 "pdb_backup": "pdb_2"},
             "set_end1":        {"src_files": "pdb_3"}}


##########################################################################
#                 Collect All Results & Output                           #
##########################################################################
class CollectResults(object):
    def __init__(self, **kwargs):
        self.steps = ["set_end", "tar_out"]

        self.recipe = \
            {"set_end":         {"command":    "set_stage_init",  # 1
                                 "options":   {"src_dir":    "",
                                               "src_files": ["pymodsim.log"],
                                               "tgt_dir":    "finalOutput"}},
                
             "tar_out":         {"command": "tar_out",  # 2
                                "options": {"src_dir": "finalOutput",
                                            "tgt":     "Model_output.tgz"}}}
            
        self.breaks = {}