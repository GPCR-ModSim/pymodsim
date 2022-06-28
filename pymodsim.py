#!/usr/bin/env python3
import argparse
import os
import shutil

import commands
import settings


class Run(object):
    def __init__(self, *args, **kwargs):
        """
        Prepare input to run PyModSim. Depending on the steps selected,
        the input will be processed accordingly.
        """
        self.own_dir = kwargs.get("own_dir")
        self.repo_dir = kwargs.get("repo_dir")
        self.sequence = kwargs.get("sequence")
        self.nstep = kwargs.get("nstep")
        self.pdb = kwargs.get("pdb") or ""
        self.Nterm = kwargs.get("Nterm") or ""
        self.Cterm = kwargs.get("Cterm") or ""
        self.loop = kwargs.get("loop") or ""
        self.replace = kwargs.get("replace") or ""

        self.cmd = commands.Commands(sequence=self.sequence,
                                     pdb=self.pdb)
        
    def clean(self):
        """
        Cleans all previously generated files
        """
        # re-write log file
        f = open("pymodsim.log", "w")
        f.close()
        
        to_unlink = []
        dirs_to_unlink = ["finalOutput"]
        
        if self.nstep == "0" or self.nstep == "1":
            unlink_files = ["ld.so.cache", "ranked_0.pdb"]
            for target in unlink_files:
                to_unlink.append(target)
            unlink_dirs = [self.sequence[:-6], "sequences"]
            for target in unlink_dirs:
                dirs_to_unlink.append(target)
        
        if self.nstep == "0" or self.nstep == "2":
            unlink_files = ["alignment.pir", "refined.B99990001.pdb", 
                            "refined.D00000001", "refined.ini", "refined.rsr",
                            "refined.sch", "refined.V99990001"]
            for target in unlink_files:
                to_unlink.append(target)
   
        if self.nstep == "0" or self.nstep == "3":
            unlink_files = ["2membranes.inp", "datapar1", "datapar2", 
                            "datasub1", "fort.4", "fort.41",
                            "homology.pdb", "refined.B99990001out.pdb", 
                            "res.lib"]
            for target in unlink_files:
                to_unlink.append(target)

        for target in to_unlink:
            if os.path.isfile(target): os.unlink(target)

        for target in dirs_to_unlink:
            if os.path.isdir(target): shutil.rmtree(target)

    def modsim(self):
        """
        Run PyModSim
        """
        #TODO: create recipes
        if self.nstep == "0":
            steps = ["Homology", "ModPrep", "Orientation"]
        if self.nstep == "1":
            steps = ["Homology"]
        if self.nstep == "2":
            steps = ["ModPrep"]
        if self.nstep == "3":
            steps = ["Orientation"]

        for step in steps:
            self.cmd.select_recipe(stage=step)
            #TODO: cmd.run_recipe
            self.cmd.run_recipe()        

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='pymodsim',
        description=' == Create prepared homology models given a sequence. == ')

    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s 0.0')

    parser.add_argument('-b',
                        dest = "own_dir",
                        help = "Working dir if different from actual dir",
                        default = os.getcwd())

    parser.add_argument('-r',
                        dest = "repo_dir",
                        help = "Path to templates of fixed files. If not \
            provided, take the value from settings.TEMPLATES_DIR.",
                        default = settings.TEMPLATES_DIR)

    parser.add_argument('-s', "--seq",
                        dest = "sequence",
                        required = True,
                        help = "Name of the sequence for which to create an \
            homology model. Use the fasta extension. (e.g. -s myseq.fasta)")

    parser.add_argument('-n', "--nstep",
                        dest = "nstep",
                        help = "PyModSim steps you wish to execute. This \
            allowes you modify the model preparation steps - see documentation \
            Options: (0) Full, (1) AlphaFold, (2) ModPrep, (3) PPM",
                        default = "0")

    parser.add_argument('-p', "--pdb",
                        dest = "pdb",
                        help = "Name of the protein to prepare. -p is only required \
            if -n = 2 or 3. Use the pdb extension. (e.g. -p myprot.pdb)",
                        default = "ranked_0.pdb")
            
    parser.add_argument('-N', "--Nterm",
                        dest = "Nterm",
                        help = "Residue number at which to cut the N-terminus. \
            Note: the chain up AND including the given residue will be removed. \
            -N is only used if -n = 2. If you wish to use the default cutoff, \
            don't specify -N. If you don't wish to cut the N-terminus: set -N = 0")            
            
    parser.add_argument('-C', "--Cterm",
                        dest = "Cterm",
                        help = "Residue number at which to cut the C-terminus. \
            Note: the chain from AND including the given residue will be removed. \
            -C is only used if -n = 2. If you wish to use the default cutoff, \
            don't specify -C. If you don't wish to cut the C-terminus: set -C = 0")             

    parser.add_argument('-l', "--loop",
                        dest = "loop",
                        help = "Residue numbers at which to cut loop(s). Define \
            the first and last residue of the loop you wish to cut (space delimited) \
            If there are multiple loops to cut, delimit the loop cuts with a '|' \
            (e.g. '101 131|230 250')")
            
    parser.add_argument('-m', "--mutate",
                        dest = "replace",
                        help = "Amino Acid sequence with which to replace the \
            cut loop. (e.g. 'AAAAA')")            

    args = parser.parse_args()

    if not (os.path.isdir(args.own_dir)):
        os.makedirs(args.own_dir)
        print("Created working dir {0}".format(args.own_dir))
    os.chdir(args.own_dir)

    run = Run(own_dir = args.own_dir,
              repo_dir = args.repo_dir,
              sequence = args.sequence,
              nstep = args.nstep,
              pdb = args.pdb,
              Nterm = args.Nterm,
              Cterm = args.Cterm,
              loop = args.loop,
              replace = args.replace)

    run.clean()
    run.modsim()
