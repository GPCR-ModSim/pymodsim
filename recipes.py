##########################################################################
#                 Homology                                               #
##########################################################################

class Homology(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "alphafold", "set_stage_init2", "set_end"]
 
        self.recipe = \
            {"set_stage_init":  {"command":    "set_stage_init",  # 1
                                 "options":   {"src_dir":    "",
                                               "src_files":  "",
                                               "tgt_dir":    "sequences"}},
                             
             "alphafold":       {"alphafold":  "alphafold",  # 2
                                 "options":   {"seq": "",
                                               "max_template_date": "2021-11-01"}},
             
             "set_stage_init2": {"command":    "set_stage_init",  # 3
                                 "options":   {"src_dir":    "",
                                               "src_files": ["ranked_0.pdb"],
                                               "tgt_dir":    "."}},
             
             
             "set_end":         {"command":    "set_stage_init",  # 4
                                 "options":   {"src_dir":    "",
                                               "src_files": ["ranked_0.pdb"],
                                               "tgt_dir":    "finalOutput"}}}
                                 
        self.breaks = \
            {"set_stage_init":  {"src_files": "sequence"},
             "alphafold":       {"seq":       "sequence"},
             "set_stage_init2": {"src_dir":   "sequence_base"}}          


##########################################################################
#                 ModPrep                                                #
##########################################################################

class ModPrep(object):
    def __init__(self, **kwargs):
        self.steps = ["plot_conf", "make_pir", "run_modeller", "set_end"]
        
        self.recipe = \
            {"plot_conf":     {"command":   "plot_conf",  # 1
                               "options":  {"pdb": "",
                                            "tgt": "plot_confidence.txt"}},
                
            "make_pir":        {"command":  "make_pir",  # 2
                                "options": {"seq": "",
                                            "pdb": "",
                                            "tgt1": "alignment.pir",
                                            "tgt2": "refinement.txt",
                                            "tgt3": "alignment.txt"}},
                         
            "run_modeller":    {"command":   "run_modeller",  # 3
                                "options":  {"alnfile": "alignment.pir",
                                             "knowns":  "",
                                             "mode":    ""}},
                        
            "set_end":         {"command":    "set_stage_init",  # 4
                                "options":   {"src_dir":    "",
                                              "src_files": ["plot_confidence.txt",
                                                            "alignment.txt",
                                                            "refinement.txt",
                                                            "refined.B99990001.pdb"],
                                              "tgt_dir":    "finalOutput"}}}                                
            
        self.breaks = \
            {"plot_conf":       {"pdb": "pdb_2"},
             "make_pir":        {"seq": "sequence",
                                 "pdb": "pdb_2"},
             "run_modeller":    {"knowns": "pdb_2",
                                 "mode": "mode"}}


##########################################################################
#                 Orientation                                            #
##########################################################################

class Orientation(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "make_inp", "ppm", "clean_pdb", 
                      "set_end1", "set_end2"]
        
        self.recipe = \
            {"set_stage_init":  {"command":  "set_stage_init",  # 1
                                 "options": {"repo_files": ["res.lib"],
                                             "tgt_dir":     "."}},
             
             "make_inp":        {"command":  "make_inp",  # 2
                                 "options": {"pdb": "",
                                             "tgt": "2membranes.inp"}},
           
             "ppm":             {"ppm":      "ppm",  # 3
                                 "options": {"inp": "2membranes.inp"}},
             
             "clean_pdb":       {"command":  "clean_pdb",  # 4
                                 "options": {"src": "",
                                             "tgt": "homology.pdb"}},
             
             "set_end1":         {"command":    "set_stage_init",  # 5
                                 "options":   {"src_dir":    "",
                                               "src_files":  "",
                                               "tgt_dir":    "finalOutput"}},
             
             "set_end2":         {"command":    "set_stage_init",  # 6
                                 "options":   {"src_dir":    "",
                                               "src_files":  ["homology.pdb"],
                                               "tgt_dir":    "finalOutput"}}}
            
        self.breaks = \
            {"make_inp":        {"pdb": "pdb_3"},
             "clean_pdb":       {"src": "pdb_3out"},
             "set_end1":        {"src_files": "pdb_3out"}}


##########################################################################
#                 Collect All Results & Output                           #
##########################################################################
class CollectResults(object):
    def __init__(self, **kwargs):
        self.steps = ["tar_out"]

        self.recipe = \
            {"tar_out":         {"command": "tar_out",  # 1
                                "options": {"src_dir": "finalOutput",
                                            "tgt":     "Model_output.tgz"}}}
            
        self.breaks = {}