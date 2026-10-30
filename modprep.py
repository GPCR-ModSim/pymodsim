import modeller
from modeller import automodel
from modeller.selection import Selection
import sys


def run_modeller(**kwargs):
    """
    Refine protein structure using MODELLER. Only the edited regions move.
    """
    default_stdout = sys.stdout
    sys.stdout = open("pymodsim.log", "a")

    try:
        if kwargs["mode"] == "run":
            modeller.log.minimal()
            env = modeller.Environ()
            env.io.atom_files_directory = ['.', '../atom_files']
            aln = modeller.Alignment(env, file=kwargs["alnfile"])

            def build(ext, randomize):
                class EditedRegionsOnly(automodel.AutoModel):
                    def select_atoms(self):
                        # Linkers and cut loop ends (+ext residues each side) move; everything
                        # else keeps the input coordinates. loops() does not report trimmed
                        # termini, so the first and last ext+1 residues always move.
                        res = self.residues
                        n = len(res)
                        ends = [res[i] for i in range(ext + 1)] + [res[n - 1 - i] for i in range(ext + 1)]
                        return Selection(self.loops(aln, minlength=1, maxlength=10000,
                                                    insertion_ext=ext, deletion_ext=ext,
                                                    include_termini=True), ends)

                a = EditedRegionsOnly(env, alnfile=kwargs["alnfile"],
                                      knowns=kwargs["knowns"], sequence='refined')
                a.starting_model = 1
                a.ending_model = 1
                a.md_level = None
                if not randomize:
                    a.rand_method = None
                a.make()
                return a.outputs[0]["failure"] is None

            # The random start can make the objective function NaN when few atoms are free
            # (2 of 1601 GPCRdb models). Retry without randomization, then with a wider selection.
            for ext, randomize in ((1, True), (1, False), (3, False)):
                if build(ext, randomize):
                    break
                print("MODELLER failed (selection +{0}, randomize={1}); retrying".format(ext, randomize))
            else:
                raise RuntimeError("MODELLER produced no model (see pymodsim.log)")
    finally:
        sys.stdout = default_stdout
