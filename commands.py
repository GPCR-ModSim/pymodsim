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

        logging.basicConfig(filename='pymodsim.log',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S',
                            level=logging.DEBUG)

        if "sequence" in kwargs.keys():
            self.set_sequence(kwargs["sequence"])
            self.set_sequence_base(kwargs["sequence"][:-6])
            
        if "pdb" in kwargs.keys():
            self.set_pdb(kwargs["pdb"])

    def set_pdb(self, value):
        """
        set_pdb: Sets the pdb object
        """
        self._pdb = value
        
    def get_pdb(self):
        return self._pdb
    pdb = property(get_pdb, set_pdb)

    def set_sequence(self, value):
        """
        set_sequence: Sets the sequence object
        """
        self._sequence = value

    def get_sequence(self):
        return self._sequence
    sequence = property(get_sequence, set_sequence)

    def set_sequence_base(self, value):
        """
        set_sequence_base: Sets the sequence object
        """
        self._sequence_base = value

    def get_sequence_base(self):
        return self._sequence_base
    sequence_base = property(get_sequence_base, set_sequence_base)

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
        tgt.write("OPC\n")                  # membrane type (DOPC bilayer)
        tgt.write("planar\n")               # flat membrane
        tgt.write("out\n")                  # N-term extracellular
        tgt.write("A\n")                    # subunits in membrane  
        
        tgt.close()

    def make_pir(self, **kwargs):
        """
        make_pir: Identify and modify low-confidence regions and create a \
        MODELLER alignment file (.pir)
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
                
        terms, loops = self.make_pir_identify(lines_pdb)  
        if len(terms) or len(loops) != 0:
            self.mode = "run"
        
        mod_seq, tmpl_seq = self.make_pir_modify(terms, loops, lines_pdb, sequence)

        tgt = open(os.path.join(self.work_dir, kwargs["tgt"]), "w") 
        
        tgt.write(f'>P1;{kwargs["pdb"]}\n')
        tgt.write(f'structure:{kwargs["pdb"]}:FIRST:@:LAST:@::::\n')

        new_line = 0
        for aa in tmpl_seq:
            tgt.write(aa)            
            new_line += 1
            if (new_line % 60) == 0:
                tgt.write("\n")

        tgt.write("\n>P1;refined\n")
        tgt.write("sequence:::::::::\n")

        new_line = 0
        for aa in mod_seq:
            tgt.write(aa)  
            new_line += 1
            if (new_line % 60) == 0:
                tgt.write("\n")

        tgt.close()

    def make_pir_identify(self, lines_pdb):
        chain = []
        low_confs = []
        loops = []
        terms = []
        
        nterm_line = lines_pdb[0].split()
        cterm_line = lines_pdb[len(lines_pdb)-3].split()
        nterm = int(nterm_line[5])
        cterm = int(cterm_line[5])

        for line in lines_pdb:
            line = line.split()
            # Only res with low model confidence (pLDDT <= 70)
            if line[0] == "ATOM" and float(line[10]) <= float(70):
                if int(line[5]) - 1 in chain:
                    if int(line[5]) not in chain:
                        chain.append(int(line[5]))
                else:
                    if chain and chain not in low_confs:
                        low_confs.append(chain)
                    if int(line[5]) not in chain:
                        chain = [int(line[5])]
                        
        for low_conf in low_confs:
            # N-term removal if low_conf longer than 5 residues
            if nterm in low_conf:
                if len(low_conf) >= 6:
                    terms.append(low_conf[:-5])
                    self.broker.dispatch("Low-confidence N-term detected from {0} to {1}".format(
                        low_conf[0], low_conf[-1]))
            #C-term removal if low_conf longer than 5 residues
            elif cterm in low_conf:
                if len(low_conf) >= 6:
                    terms.append(low_conf[5:])
                    self.broker.dispatch("Low-confidence C-term detected from {0} to {1}".format(
                        low_conf[0], low_conf[-1]))
            #Loop removal if low_conf longer than 10 residues
            else: 
                if len(low_conf) >= 11:
                    loops.append(low_conf)
                    self.broker.dispatch("Low-confidence loop detected from {0} to {1}".format(
                        low_conf[0], low_conf[-1]))
    
        loops.reverse()

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
                line = line.split()
                if line[0] == "ATOM":
                    if int(line[5]) == loop[0] and line[2] == "C":
                        start = [float(line[6]), float(line[7]), float(line[8])]
                    if int(line[5]) == loop[-1] and line[2] == "N":
                        end = [float(line[6]), float(line[7]), float(line[8])]

            # AA countour length is between 3.4 and 4.0 Å. To allow for 
            # flexibility in the loop, 1 AA per 2.5 Å is used. 
            # TODO: adjust AA per Å to see effects
            aa_dist = 2.0
            
            x = abs(start[0] - end[0])
            y = abs(start[1] - end[1])
            z = abs(start[2] - end[2])
            dist = math.sqrt(x*x + y*y + z*z)
            loop_fill = math.ceil(dist / aa_dist)     
            
            if loop_fill % 2 == 0:
                keep_res = int(loop_fill / 2 - 1)
                ala_count = 2
            else: 
                keep_res = int((loop_fill - 1) / 2 - 1) 
                ala_count = 3
            
            # Leaving the first and the last AA out of the loop (for better results in 
            # MODELLER)
            loop_seq = sequence[loop[0]:loop[-2]]
            loop_mod = loop_seq[:keep_res-1] + "A"*ala_count + loop_seq[-1*(keep_res-1):]
            gaps = loop_mod + len(loop_seq) * "-"
            
            mod_seq = mod_seq.replace(loop_seq, gaps)
            tmpl_seq = sequence[:loop[0]] + len(loop_mod) * "-" + sequence[loop[0]:]
            
        mod_seq += "*"
        tmpl_seq += "*"

        return mod_seq, tmpl_seq
       
    def run_recipe(self):
        #TODO: currently modified
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
        self.recipe = getattr(recipes, stage)(sequence=self.sequence)
         
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
        command.extend(["--max_template_date=2021-11-01"])
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