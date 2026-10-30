PyModSim Version 1.0
================================================================================

PyModSim is a standalone *Python* package to create protein homology models from 
a sequence using **AlphaFold** (or **AlphaFold MultiState**, to model a GPCR in 
its active or inactive state), refine protein models using **MODELLER**, and 
standardize the orientation of protein models using **PPM**.

**PyModSim** is hosted in GitHub at:

<https://github.com/rlvandenbroek/pymodsim>

You can download any version of **PyModSim** by cloning the repository to your 
local machine using git.  

You will need to create a  free personal account at github and send
and  e-mail  to:  [r.l.van.den.broek@lacdr.leidenuniv.nl](r.l.van.den.broek@lacdr.leidenuniv.nl) 
requesting access to the code. After request processing from us you will be
given access to the free repository.  

To install **PyModSim** follow these steps:  

0.  Install/prepare all **PyModSim** dependencies:
	
        Python 3.X (extensively tested on Python 3.10)
        Python modules:
            - modeller 10.2 (https://salilab.org/modeller/release.html)
            - biopython 1.81 (https://biopython.org/wiki/Download)
        AlphaFold 2.0 - (https://github.com/deepmind/alphafold)
        AlphaFold MultiState - only for -n 10 and 11 (https://github.com/huhlim/alphafold-multistate)
        PPM 3.0 - (https://console.cloud.google.com/storage/browser/opm-assets/ppm3_code)

1.  Clone **PyModSim** for Python 3.X with the *modeller* module:  

        git clone https://username@github.com/GPCR-ModSim/pymodsim.git

    Make sure to change *username* to the one you have created at
    github.  

2.  The previous command will create a *pymodsim* directory. Now you
    have to tell your operating system how to find that folder. You
    achieve this by declaring the location of the directory in a .bashrc
    file .cshrc or .zshrc file in your home folder. An example of what you will
    have to include in your .bashrc file follows:

        export PYMODSIM=/home/username/software/pymodsim
        export PATH=$PYMODSIM:$PATH

    or if your shell is csh then in your .cshrc file you can add:

        setenv PYMODSIM /home/username/software/pymodsim
        set path = ($path $PYMODSIM)

    Notice that I have cloned *pymodsim* in the software folder in my
    home folder, you will have to adapt this to wherever it is that you
    downloaded your *pymodsim* to.

    After including the route to your *pymodsim* directory in your
    .bashrc file make sure to issue the command:

        source .bashrc

    or open a new terminal.

    To check if you have defined the route to the *pymodsim* directory
    correctly try to run the main program called pymodsim in a terminal:

        pymodsim -h

    You should obtain the following help output:
	
        usage: pymodsim [-h] [-v] [-n NSTEP] [-s SEQUENCE] [-p PDB] [-N NTERM] [-C CTERM] [-l LOOP]
                [-f LOOP_FILL] [-t TOPOLOGY] [-c CHAIN] [-a {active,inactive}]
	
        == Create prepared homology models given a sequence. ==

        options:
          -h, --help            show this help message and exit
          -v, --version         show program's version number and exit
          -n NSTEP, --nstep NSTEP
                                PyModSim steps you wish to execute. This allowes you modify the
                                model preparation steps - see documentation Options: (0) Full, (1)
                                Homology|AlphaFold, (2) ModPrep|MODELLER, (3) Orientation|PPM, (23)
                                ModPrep+Orientation, (10) Full with StateHomology|AlphaFold
                                MultiState, and (11) StateHomology|AlphaFold MultiState. -n 10 and
                                11 require -s and -a
          -s SEQUENCE, --seq SEQUENCE
                                Name of the sequence for which to create an homology model. -s is
                                only required if -n = 0, 1, 2, 23, 10 or 11. Use the fasta
                                extension. (example: -s myseq.fasta)
          -p PDB, --pdb PDB     Name of the protein to process. -p is only required if -n = 2 or 3.
                                Use the pdb extension. (example: -p myprot.pdb)
          -N NTERM, --Nterm NTERM
                                Residue number at which to cut the N-terminus. Note: the chain up
                                AND including the given residue will be removed. -N is only used if
                                -n = 2. If you wish to use the default cutoff, don't specify -N. If
                                you wish not to cut the N-term: set -N = 0
          -C CTERM, --Cterm CTERM
                                Residue number at which to cut the C-terminus. Note: the chain from
                                AND including the given residue will be removed. -C is only used if
                                -n = 2. If you wish to use the default cutoff, don't specify -C. If
                                you wish not to cut the C-term: set -C = 0
          -l LOOP, --loop LOOP  Residue numbers at which to cut loop(s). Define the first and last
                                residue of the loop you wish to cut ('-' delimited) If there are
                                multiple loops to cut, delimit the loop cuts with a ',' (example: -l
                                101-131,230-250). If you do not with to cut any loops: set -l = 0
          -f LOOP_FILL, --loop_fill LOOP_FILL
                                Amount of Å per AA to fill cut loops. The total distance is
                                calculated from the coordinates of the remaining residues. The AA
                                contour length is 3.4-4.0 Å, To allow for flexibility in the loop,
                                2.0 Å/AA (default) is suggested. (example: -f 2.0)
          -t TOPOLOGY, --topology TOPOLOGY
                                Indicate the topology of the N-term within the protein structure.
                                'out': extracellular N-term (default), 'in': intracellular N-term.
          -c CHAIN, --chain CHAIN
                                Only use if -n = 3 (i.e. only PPM). If more than 1 chain, add a
                                comma-seperated list of the chain identifiers. (example: -c A,B,C)
          -a {active,inactive}, --state {active,inactive}
                                Receptor state to model with AlphaFold MultiState (GPCR100 Active
                                or Inactive templates). -a is required if -n = 10 or 11. With -n 23
                                it is optional and only adds a copy of the final model named
                                <state>_homology.pdb. (example: -a active)

3.  Updates are very easy thanks to the git versioning system. Once
    **PyModSim** has been downloaded (cloned) into its own *pymodsim* folder 
    you just have to move to it and pull the newest changes:

        cd /home/username/software/pymodsim
        git pull   

5.  To make sure that your AlphaFold and PPM installation is understood by
    **PyModSim** you will need to specify the path to where the sofware is
    installed in your system. To do this you will need to edit the
    settings.py file with any text editor (“vi” and “emacs” are common
    options in the unix environment). Make sure that only one line is
    uncommented, looking like: PPM_PATH = /apps/PPM_3.0 Provided that in 
    your case PPM is installed in /apps. The program
    will prepend this line to the binaries names, so calling
    “/apps/PPM_3.0/immers should point to that binary.  

    For -n 10 and 11, also set AFMS_SCRIPT in settings.py: the launcher that
    runs AlphaFold MultiState with its own environment and databases. PyModSim
    calls it as:

        AFMS_SCRIPT <fasta> <output_dir> GPCR <active|inactive>

    It must exit with 0 and write <output_dir>/<active|inactive>_model.pdb.
    AFMS_UNSET_ENV lists environment variables that are removed before calling
    it, for example the MODELLER paths, which are not meant for the AlphaFold
    MultiState environment, whose script, at the time of writing, was in the
    following path:

        AFMS_SCRIPT = "/path/to/alphafold_run.sh"
        AFMS_UNSET_ENV = ["PYTHONPATH", "LD_LIBRARY_PATH"]

    Modes 0, 1, 2, 3 and 23 do not need these settings.

### Modelling a GPCR in a given state (-n 10 and 11)

The StateHomology step uses AlphaFold MultiState instead of AlphaFold. Its
templates are annotated by activation state (GPCR100 Active or Inactive), so a
GPCR can be modelled in the active or inactive state. Choose the state with -a:

        pymodsim -n 10 -s myseq.fasta -a active      # StateHomology + ModPrep + Orientation
        pymodsim -n 11 -s myseq.fasta -a inactive    # StateHomology only (raw model)

- AlphaFold MultiState runs in myseq_afms/, and its output is written to
  pymodsim.log while it runs. Its model, myseq_afms/<state>_model.pdb, is
  copied to ranked_0.pdb, so ModPrep and Orientation process it exactly like
  an AlphaFold model.
- Model_output.tgz contains the usual files plus <state>_homology.pdb
  (-n 10, prepared and oriented) or <state>_model.pdb (-n 11, raw model).
- If AlphaFold MultiState fails, PyModSim stops with exit code 1 and writes no
  Model_output.tgz. Invalid arguments stop it with exit code 2 before any file
  is deleted.
- A prediction can take many hours. If the preparation fails afterwards,
  rerun only ModPrep + Orientation on the existing model, instead of repeating
  -n 10, which deletes myseq_afms/:

        pymodsim -n 23 -s myseq.fasta -p ranked_0.pdb -c A -a active

### Auxiliary Modules

- **recipes.py**.   Applies step by step instructions for carrying a 
  modeling step (Homology, StateHomology, ModPrep, Orientation, CollectResults).
- **broker.py**.   Proxy for printing messages
- **settings.py**.   This modules sets up the main environment variables needed
  to run the calculation, for example, the path to the AlphaFold and PPM binaries
  and to the AlphaFold MultiState launcher (AFMS_SCRIPT).

### Execution Modules

- **commands.py**. Defines the Commands and Wrapper objects. Commands will
  load the objects, recipes, and run them. Wrapper is a proxy for the 
  commands. When a recipe entry is sent to it this returns the command to be run.
- **modprep.py**. Runs commands using the *modeller* module.

### Executable

- **pymodsim** The main program to call which sends the run to a cluster.
