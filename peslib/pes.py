import numpy as np
from ase import Atoms
from ase.units import Bohr, Angstrom, Hartree, eV
from peslibf import o4_singlet
from ase.calculators.calculator import Calculator, all_changes


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
        r = atoms.get_positions()
        x = r[:, 0] / Bohr * Angstrom
        y = r[:, 1] / Bohr * Angstrom
        z = r[:, 2] / Bohr * Angstrom
        dx = np.zeros(4)
        dy = np.zeros(4)
        dz = np.zeros(4)
        e = np.zeros(2)
        o4_singlet.pot(x, y, z, e, dx, dy, dz)
        f = np.zeros((4, 3))
        f[:, 0] = -dx
        f[:, 1] = -dy
        f[:, 2] = -dz
        e = e[0] * Hartree / eV
        self.results = {'energy': e, 'forces': f}


if __name__ == '__main__':
    atoms = Atoms('O4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    atoms.calc = O4SingletPES()
    atoms.get_potential_energy()
