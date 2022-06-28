class Homology(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "alphafold", "set_stage_init2", "set_end"]
 
        self.recipe = \
            {"set_stage_init":  {"command":    "set_stage_init",  # 1
                                 "options":   {"src_dir":    "",
                                               "src_files":  "",
                                               "tgt_dir":    "sequences"}},
                             
             "alphafold":       {"alphafold":  "alphafold",  # 2
                                 "options":   {"seq": ""}},
             
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


class ModPrep(object):
    def __init__(self, **kwargs):
        self.steps = ["make_pir", "run_modeller", "set_end"]
        
        self.recipe = \
            {"make_pir":       {"command":  "make_pir",  # 1
                                "options": {"seq": "",
                                            "pdb": "",
                                            "tgt": "alignment.pir"}},
                         
            "run_modeller":    {"command":   "run_modeller",  # 2
                                "options":  {"alnfile": "alignment.pir",
                                             "knowns":  "ranked_0.pdb",
                                             "mode": ""}},
                        
            "set_end":         {"command":    "set_stage_init",  # 3
                                "options":   {"src_dir":    "",
                                              "src_files": ["refined.B99990001.pdb"],
                                              "tgt_dir":    "finalOutput"}}}                                
            
        self.breaks = \
            {"make_pir":        {"seq": "sequence",
                                 "pdb": "pdb"},
             "run_modeller":    {"mode": "mode"}}


class Orientation(object):
    def __init__(self, **kwargs):
        self.steps = ["set_stage_init", "make_inp", "ppm", "clean_pdb", 
                      "set_end"]
        
        self.recipe = \
            {"set_stage_init":  {"command":  "set_stage_init",  # 1
                                 "options": {"repo_files": ["res.lib"],
                                             "tgt_dir":     "."}},
             
             "make_inp":        {"command":  "make_inp",  # 2
                                 "options": {"pdb": "refined.B99990001.pdb",
                                             "tgt": "2membranes.inp"}},
           
             "ppm":             {"ppm":      "ppm",  # 3
                                 "options": {"inp": "2membranes.inp"}},
             
             "clean_pdb":       {"command":  "clean_pdb",  # 4
                                 "options": {"src": "refined.B99990001out.pdb",
                                             "tgt": "homology.pdb"}},
             
             "set_end":         {"command":    "set_stage_init",  # 5
                                 "options":   {"src_dir":    "",
                                               "src_files": ["refined.B99990001out.pdb",
                                                             "homology.pdb"],
                                               "tgt_dir":    "finalOutput"}}}
            
        self.breaks = {}