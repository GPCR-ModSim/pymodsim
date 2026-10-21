import modeller
from modeller import automodel
from modeller.selection import Selection
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
        aln = modeller.Alignment(env, file=kwargs["alnfile"])

        class EditedRegionsOnly(automodel.AutoModel):
            def select_atoms(self):
                # Linkers and cut loop ends (+1 residue each side) move; everything else keeps
                # the input coordinates. loops() does not report trimmed termini, so the first
                # and last two residues always move.
                res = self.residues
                n = len(res)
                return Selection(self.loops(aln, minlength=1, maxlength=10000,
                                            insertion_ext=1, deletion_ext=1,
                                            include_termini=True),
                                 [res[0], res[1], res[n - 2], res[n - 1]])

        a = EditedRegionsOnly(env, alnfile=kwargs["alnfile"],
                              knowns=kwargs["knowns"], sequence='refined')
        a.starting_model = 1
        a.ending_model = 1
        a.md_level = None
        a.make()

    sys.stdout = default_stdout
