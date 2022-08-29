import functools
import logging
import math
import os
import shutil
import subprocess

import broker
import modprep
import recipes
import settings


class Commands(object):
    def __init__(self, *args, **kwargs):
        self.broker = broker.Printing()
        self.wrapper = Wrapper()        
        self.work_dir = os.getcwd()
        self.repo_dir = self.wrapper.repo_dir
        self.mode = "skip"        

        if kwargs["sequence"]:
            self.sequence = kwargs["sequence"]
            self.sequence_base = kwargs["sequence"][:-6]
        if kwargs["pdb"]:
            self.pdb_2 = kwargs["pdb"]
        else:
            self.pdb_2 = "ranked_0.pdb"            
        if kwargs["pdb"] and kwargs["nstep"] != "23":
            self.pdb_3 = kwargs["pdb"]
            self.pdb_3out = kwargs["pdb"][:-4] + "out.pdb"
        else:
            self.pdb_3 = "refined.B99990001.pdb"
            self.pdb_3out = "refined.B99990001out.pdb"
        self.Nterm = kwargs["Nterm"]
        self.Cterm = kwargs["Cterm"]
        self.loop = kwargs["loop"]
        self.loop_fill = kwargs["loop_fill"]
        self.topology = kwargs["topology"]
        self.chain = kwargs["chain"]
        
        logging.basicConfig(filename='pymodsim.log',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S',
                            level=logging.DEBUG)

    def clean_pdb(self, **kwargs):
        """
        clean_pdb: Remove membrane from pdb
        """
        src = open(os.path.join(self.work_dir, kwargs["src"]), "r")
        lines_src = src.readlines()
        src.close()
        
        tgt = open(os.path.join(self.work_dir, kwargs["tgt"]), "w") 

        for line in lines_src:
            if line[0:6] != "HETATM":
                tgt.write(line)
                
        tgt.close()

    def make_inp(self, **kwargs):
        """
        make_inp: Create a PPM input file (.inp)
        """
        tgt = open(os.path.join(self.work_dir, kwargs["tgt"]), "w")
        
        tgt.write("2\n")                    # PPM mode
        tgt.write("no\n")                   # no heteroatoms
        tgt.write(kwargs["pdb"] + "\n")     # PDB file name
        tgt.write("1\n")                    # number of membranes
        tgt.write("PMm\n")                  # membrane type (Plasma membrane (mammalian))
        tgt.write("planar\n")               # flat membrane
        tgt.write(self.topology + "\n")     # N-term topology
        tgt.write(self.chain + "\n")        # subunits in membrane  

        tgt.close()

    def make_pir(self, **kwargs):
        """
        make_pir: Identify and modify low-confidence regions and create a MODELLER alignment file (.pir)
        """
        pdb = open(os.path.join(self.work_dir, kwargs["pdb"]), "r")
        seq = open(os.path.join(self.work_dir, kwargs["seq"]), "r")
        lines_pdb = pdb.readlines()
        lines_seq = seq.readlines()
        pdb.close()
        seq.close()
        
        sequence = ""
        for line in lines_seq:
            if line[0] != ">":
                sequence += line.replace("\n", "")
                
        terms, loops = self.make_pir_identify(lines_pdb, kwargs["tgt2"])  
        if len(terms) or len(loops) != 0:
            self.mode = "run"
        
        mod_seq, tmpl_seq = self.make_pir_modify(terms, loops, lines_pdb, sequence)

        tgt1 = open(os.path.join(self.work_dir, kwargs["tgt1"]), "w") 
        tgt3 = open(os.path.join(self.work_dir, kwargs["tgt3"]), "w")
        
        tgt1.write(f'>P1;{kwargs["pdb"]}\n')
        tgt1.write(f'structure:{kwargs["pdb"]}:FIRST:@:LAST:@::::\n')

        new_line = 0
        for aa in tmpl_seq:
            tgt1.write(aa)            
            new_line += 1
            if (new_line % 60) == 0:
                tgt1.write("\n")

        tgt1.write("\n>P1;refined\n")
        tgt1.write("sequence:::::::::\n")

        new_line = 0
        for aa in mod_seq:
            tgt1.write(aa)  
            new_line += 1
            if (new_line % 60) == 0:
                tgt1.write("\n")

        tgt3.write("tmpl_ID\ttmpl_AA\tref_ID\tref_AA\n")
        tmpl_count = ref_count = 0
        for i in range(len(tmpl_seq)-1):
            tmpl_ID = ref_ID = "-"
            if tmpl_seq[i] != "-":
                tmpl_count += 1
                tmpl_ID = tmpl_count
            if mod_seq[i] != "-":
                ref_count += 1
                ref_ID = ref_count
            tgt3.write(str(tmpl_seq[i]) + "\t" + str(tmpl_ID) + "\t" + \
                       str(mod_seq[i]) + "\t" + str(ref_ID) + "\n")

        tgt1.close()
        tgt3.close()

    def make_pir_identify(self, lines_pdb, tgt):
        tgt2 = open(os.path.join(self.work_dir, tgt), "w")
        
        chain = []
        low_confs = []
        loops = []
        terms = []
        
        nterm_line = lines_pdb[0]
        cterm_line = lines_pdb[len(lines_pdb)-3]
        nterm = int("".join(nterm_line[22:26]))
        cterm = int("".join(cterm_line[22:26]))

        for line in lines_pdb:
            # Only res with low model confidence (pLDDT <= 70)
            if "".join(line[:4]) == "ATOM" and float("".join(line[60:66])) <= float(70):
                if int("".join(line[22:26])) - 1 in chain:
                    if int("".join(line[22:26])) not in chain:
                        chain.append(int("".join(line[22:26])))
                else:
                    if chain and chain not in low_confs:
                        low_confs.append(chain)
                    if int("".join(line[22:26])) not in chain:
                        chain = [int("".join(line[22:26]))]
        
        if self.Nterm:
            if self.Nterm != "0":
                terms.append(range(nterm, int(self.Nterm)))
                tgt2.write("Custom N-term removed from {0} to {1}\n".format(
                    nterm, self.Nterm))
                self.broker.dispatch("Custom N-term removed from {0} to {1}".format(
                    nterm, self.Nterm))                    
        
        if self.Cterm:
            if self.Cterm != "0":
                terms.append(range(int(self.Cterm), cterm+1))
                tgt2.write("Custom C-term removed from {0} to {1}\n".format(
                    self.Cterm, cterm))                 
                self.broker.dispatch("Custom C-term removed from {0} to {1}".format(
                    self.Cterm, cterm)) 
        
        if self.loop:
            if self.loop != "0":
                cust_loops = self.loop.split(",") 
                for loop_ends in cust_loops:
                    bound = loop_ends.split("-")
                    cust_loop = range(int(bound[0]), int(bound[1])+1)
                    loops.append(cust_loop)
                    tgt2.write("Custom loop removed from {0} to {1}\n".format(
                        bound[0], bound[1]))
                    self.broker.dispatch("Custom loop removed from {0} to {1}".format(
                        bound[0], bound[1]))
            
        for low_conf in low_confs:
            # N-term removal if low_conf longer than 5 residues
            incl_loop = True
            
            if not self.Nterm:
                if nterm in low_conf:
                    if len(low_conf) >= 6:
                        terms.append(low_conf[:-5])
                        tgt2.write("Low-confidence N-term detected from {0} to {1}\n".format(
                            low_conf[0], low_conf[-1]))                      
                        self.broker.dispatch("Low-confidence N-term detected from {0} to {1}".format(
                            low_conf[0], low_conf[-1]))
            
            #C-term removal if low_conf longer than 5 residues
            if not self.Cterm:
                if cterm in low_conf:
                    if len(low_conf) >= 6:
                        terms.append(low_conf[5:])
                        tgt2.write("Low-confidence C-term detected from {0} to {1}\n".format(
                            low_conf[0], low_conf[-1]))
                        self.broker.dispatch("Low-confidence C-term detected from {0} to {1}".format(
                            low_conf[0], low_conf[-1]))
               
            if not self.loop:
                if self.Nterm and self.Nterm != "0":
                    for i in list(range(nterm, int(self.Nterm))):
                        if i in low_conf:
                            incl_loop = False
                else:
                    if nterm in low_conf:
                        incl_loop = False
                
                if self.Cterm and self.Cterm != "0":
                    for j in list(range(int(self.Cterm), cterm)):
                        if j in low_conf:
                            incl_loop = False
                else:
                    if cterm in low_conf:
                        incl_loop = False
                
                if incl_loop == True:
                    if len(low_conf) >= 11:
                        loops.append(low_conf)
                        tgt2.write("Low-confidence loop detected from {0} to {1}\n".format(
                            low_conf[0], low_conf[-1]))
                        self.broker.dispatch("Low-confidence loop detected from {0} to {1}".format(
                            low_conf[0], low_conf[-1]))
                    
                    
        loops.reverse()
        
        tgt2.close()

        return terms, loops
    
    def make_pir_modify(self, terms, loops, lines_pdb, sequence):
        mod_seq = sequence
        tmpl_seq = sequence
        
        for term in terms:
            term_seq = sequence[term[0]-1:term[-1]]
            gaps =  len(term_seq) * "-"
            
            mod_seq = mod_seq.replace(term_seq, gaps)
        
        for loop in loops:                       
            for line in lines_pdb:
                if line[:4] == "ATOM":
                    if int(line[22:26]) == loop[0] and "".join(line[11:17]).strip() == "C":
                        start = [float(line[31:38]), float(line[39:46]), float(line[47:54])]
                    if int(line[22:26]) == loop[-1] and "".join(line[11:17]).strip() == "N":
                        end = [float(line[31:38]), float(line[39:46]), float(line[47:54])]
        
            aa_dist = float(self.loop_fill)           
            x = abs(start[0] - end[0])
            y = abs(start[1] - end[1])
            z = abs(start[2] - end[2])
            dist = math.sqrt(x*x + y*y + z*z)
            num_aa = math.ceil(dist / aa_dist)     
            
            if num_aa % 2 == 0:
                keep_res = int(num_aa / 2 - 1)
                ala_count = 2
            else: 
                keep_res = int((num_aa - 1) / 2 - 1) 
                ala_count = 3
            
            # Leaving the first and the last AA out of the loop (for better results in 
            # MODELLER)
            loop_seq = sequence[loop[0]:loop[-2]]
            loop_mod = loop_seq[:keep_res-1] + "A"*ala_count + loop_seq[-1*(keep_res-1):]
            gaps = loop_mod + len(loop_seq) * "-"
            
            mod_seq = mod_seq.replace(loop_seq, gaps)
            tmpl_seq = tmpl_seq[:loop[0]] + len(loop_mod) * "-" + tmpl_seq[loop[0]:]
            
        mod_seq += "*"
        tmpl_seq += "*"

        return mod_seq, tmpl_seq
      
    def plot_conf(self, **kwargs):
        """
        plot_conf: Create an AF confidence plot (res_ID vs pLDDT)
        """
        pdb = open(kwargs["pdb"], "r")
        lines = pdb.readlines()
        
        res_list = []
        conf_list = []

        for line in lines:
            try: 
                pLDDT = line[60:66]
                res = line[22:26]
                if res not in res_list:
                    res_list.append(res)
                    conf_list.append(float(pLDDT))

            except: IndexError

        conf_vh = conf_h = conf_l = conf_vl = 0

        for pLDDT in conf_list:
            if pLDDT >= 90:
                conf_vh += 1
            elif pLDDT >= 70:
                conf_h += 1
            elif pLDDT >= 50:
                conf_l += 1
            else:
                conf_vl += 1

        avg_pLDDT = round(sum(conf_list) / len(conf_list), 2)     

        per_vh = round(conf_vh / len(conf_list) * 100, 1)
        per_h = round(conf_h / len(conf_list) * 100, 1)
        per_l = round(conf_l / len(conf_list) * 100, 1)
        per_vl = round(conf_vl / len(conf_list) * 100, 1)

        tgt = open(kwargs["tgt"], "w")

        tgt.write("Confidence Report: " + str(kwargs["pdb"]) + "\n\n" + \
                  "Average pLDDT: " + str(avg_pLDDT) + "\n\n" + \
                  "Model Confidence:\n" + \
                  "Very high (pLDDT > 90): " + str(per_vh) + "%\n" + \
                  "High (90 > pLDDT > 70): " + str(per_h) + "%\n" + \
                  "Low (70 > pLDDT > 50):  " + str(per_l) + "%\n" + \
                  "Very low (pLDDT < 50):  " + str(per_vl) + "%\n\n" + \
                  "res_ID\tpLDDT\n")
            
        for i in range(len(conf_list)):
            tgt.write(str(res_list[i]) + "\t" + str(conf_list[i]) + "\n")          
            
        pdb.close()
        tgt.close()

      
    def run_recipe(self):
        """
        run_recipe: Run selected recipes 
        """        
        for n, command_name in enumerate(self.recipe.steps):
            command = self.recipe.recipe[command_name]
            
            if command_name in self.recipe.breaks.keys():
                command["options"] = self.set_options(command["options"],
                                     self.recipe.breaks[command_name])            
            
            self.broker.dispatch("{0} Step ({1}/{2}): {3}.".format(
                self.recipe.__class__.__name__,
                n + 1, len(self.recipe.steps),
                command_name))
            
            if "alphafold" in command:
                out, err = self.wrapper.run_command(prgm="alphafold", cmd=command)
                command, stdin = self.wrapper.generate_command(prgm="alphafold", cmd=command)                
                logging.debug(" ".join(command))
                logging.debug(err.decode().strip('\n'))
                logging.debug(out.decode().strip('\n'))
                
            elif "ppm" in command:
                out, err = self.wrapper.run_command(prgm="ppm", cmd=command)
                command, stdin = self.wrapper.generate_command(prgm="ppm", cmd=command)
                logging.debug(command + " < " + stdin.name)
                logging.debug(err.decode().strip('\n'))
                logging.debug(out.decode().strip('\n'))
            
            else:
                # ...or run a local function
                logging.debug(command)
                try: 
                    f = getattr(self, command["command"])
                except AttributeError:
                    f = getattr(modprep, command["command"])
                
                logging.debug("FUNCTION: " + str(f.__doc__).strip())                
                if ("options") in command:
                    f(**command["options"])
                else:
                    f()    

    def select_recipe(self, stage):
        """
        select_recipe: Select the recipes for each step
        """
        self.recipe = getattr(recipes, stage)()
         
    def set_options(self, options, breaks):
        """
        set_options: Set break options from recipe
        """
        for option, value in breaks.items():
            # This is a hack to get the attribute recursively,
            # feeding getattr with dot-splitted string thanks to reduce
            # Here we charge some commands with options calculated
            new_option = functools.reduce(getattr,
                                value.split("."),
                                self)
            options[option] = new_option

        return options

    def set_stage_init(self, **kwargs):
        """
        set_stage_init: Copy a set of files from source to target dir
        """
        if not os.path.isdir(kwargs["tgt_dir"]): os.mkdir(kwargs["tgt_dir"])

        if "src_files" in kwargs.keys():
            if not isinstance(kwargs["src_files"], list):
                kwargs["src_files"] = kwargs["src_files"].split(" ")
            
            for src_file in kwargs["src_files"]:
                if (os.path.isfile(os.path.join(kwargs["src_dir"], src_file))):
                    shutil.copy(os.path.join(kwargs["src_dir"], src_file),
                                os.path.join(kwargs["tgt_dir"],
                                         os.path.split(src_file)[1]))

        if "repo_files" in kwargs.keys():
            for repo_file in kwargs["repo_files"]:
                shutil.copy(os.path.join(self.repo_dir, repo_file),
                            os.path.join(kwargs["tgt_dir"], repo_file))


class Wrapper(object):
    def __init__(self, *args, **kwargs):
        self.work_dir = os.getcwd()
        self.repo_dir = settings.TEMPLATES_DIR
        self.cmd_type = "default"

        #alphafold variables
        self.ALPHAFOLD_SIF = settings.ALPHAFOLD_SIF
        self.ALPHAFOLD_DATA_PATH = settings.ALPHAFOLD_DATA_PATH
        self.ALPHAFOLD_MODELS = settings.ALPHAFOLD_MODELS
        self.BFD_DATABASE_PATH = settings.BFD_DATABASE_PATH
        self.MGNIFY_DATABASE_PATH = settings.MGNIFY_DATABASE_PATH
        self.PDB70_DATABASE_PATH = settings.PDB70_DATABASE_PATH
        self.TEMPLATE_MMCIF_DIR = settings.TEMPLATE_MMCIF_DIR
        self.OBSOLETE_PDBS_PATH = settings.OBSOLETE_PDBS_PATH
        self.UNICLUST30_DATABASE_PATH = settings.UNICLUST30_DATABASE_PATH
        self.UNIREF90_DATABASE_PATH = settings.UNIREF90_DATABASE_PATH
        
        #PPM variables:
        self.PPM_PATH = settings.PPM_PATH

    def generate_command(self, prgm, cmd):
        """
        generate_command: Receive some variables in kwargs, generate
        the appropriate command to be run. Return a set in the form of
        a string "command -with flags"
        """
        stdin = ""
        
        if prgm == "alphafold":
            command = self._mode_alphafold(cmd=cmd)
            
        if prgm == "ppm":
            command, stdin = self._mode_ppm(cmd=cmd)
            self.cmd_type = "input"
        
        return command, stdin

    def _mode_alphafold(self, cmd, **kwargs):
        command = ["singularity", "run"]
        command.extend(["--nv"])
        command.extend(["-B", self.ALPHAFOLD_DATA_PATH + ":/data"])
        command.extend(["-B", self.ALPHAFOLD_MODELS])
        command.extend(["-B", self.work_dir + ":/output"])
        command.extend(["-B", ".:/etc"])
        command.extend(["--pwd", "$TMPDIR", self.ALPHAFOLD_SIF])
        command.extend(["--fasta_paths=/output/sequences/" + cmd["options"]["seq"]])
        command.extend(["--data_dir=/data"])
        command.extend(["--output_dir=/output"])
        command.extend(["--bfd_database_path=" + self.BFD_DATABASE_PATH])
        command.extend(["--mgnify_database_path=" + self.MGNIFY_DATABASE_PATH])
        command.extend(["--pdb70_database_path=" + self.PDB70_DATABASE_PATH])
        command.extend(["--template_mmcif_dir=" + self.TEMPLATE_MMCIF_DIR])
        command.extend(["--obsolete_pdbs_path=" + self.OBSOLETE_PDBS_PATH])
        command.extend(["--uniclust30_database_path=" + self.UNICLUST30_DATABASE_PATH])
        command.extend(["--uniref90_database_path=" + self.UNIREF90_DATABASE_PATH]) 
        command.extend(["--max_template_date=" + cmd["options"]["max_template_date"]])
        command.extend(["--model_preset=monomer"])
        command.extend(["--db_preset=full_dbs"])
        
        return command
    
    def _mode_ppm(self, cmd, **kwargs):       
        command = os.path.join(self.PPM_PATH, "immers")
        stdin = open(os.path.join(self.work_dir, cmd["options"]["inp"]), "r")
        
        return command, stdin
    
    def run_command(self, prgm, cmd):
        command, stdin = self.generate_command(prgm=prgm, cmd=cmd)

        if self.cmd_type == "input":
            p = subprocess.Popen(command,
                                 stdin=stdin,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen(command,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)       

        out, errs = p.communicate()

        return out, errs
    