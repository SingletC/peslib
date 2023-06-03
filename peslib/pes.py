import numpy as np
from ase import Atoms
from ase.units import Bohr, Angstrom, Hartree, eV
from peslibf import o4_singlet
from ase.calculators.calculator import Calculator, all_changes

ang2bohr = Angstrom / Bohr
hatree2ev = Hartree/eV
hatree_bohr2ev_ang = hatree2ev * ang2bohr
class O4SingletPES(Calculator):
    """
    O4_singlet
    Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,D. G. Truhlar J.Chem.Phys.
    """
    implemented_properties = [
        "energy",
        "forces", ]

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            self.atoms = atoms.copy()
        r = atoms.get_positions() * ang2bohr
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        dx = np.zeros(4)
        dy = np.zeros(4)
        dz = np.zeros(4)
        e = np.zeros(2)
        o4_singlet.pot(x, y, z, e, dx, dy, dz)
        f = np.zeros((4, 3))
        f[:, 0] = -dx * hatree_bohr2ev_ang
        f[:, 1] = -dy * hatree_bohr2ev_ang
        f[:, 2] = -dz * hatree_bohr2ev_ang
        e = e[0] * hatree2ev
        self.results = {'energy': e, 'forces': f}


if __name__ == '__main__':
    atoms = Atoms('O4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    atoms.calc = O4SingletPES()
    atoms.get_potential_energy()
