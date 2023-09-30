from typing import Callable, Tuple, List, Optional

import numpy as np
from ase import Atoms
from ase.units import Bohr, Angstrom, Hartree, eV
from peslibf import o4_singlet, n4_singlet, o4_triplet, n2o2_triplet
from ase.calculators.calculator import Calculator, all_changes

ang2bohr = Angstrom / Bohr
hatree2ev = Hartree / eV
hatree_bohr2ev_ang = hatree2ev * ang2bohr


class BasePES(Calculator):
    __pes__func__: Callable[
        [np.ndarray, np.ndarray, np.ndarray], Tuple[float, np.ndarray, np.ndarray, np.ndarray]] = None
    __atomic_numbers__: Optional[List[int]] = None

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            self.atoms = atoms.copy()
        if not (atoms.get_atomic_numbers() == self.__atomic_numbers__).all() and self.__atomic_numbers__:
            raise ValueError(
                f'order of Atomic numbers of atoms should be {self.__atomic_numbers__}')
        r = atoms.get_positions() * ang2bohr
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        e, dx, dy, dz = self.__pes__func__(x, y, z)
        f = np.empty((4, 3))
        f[:, 0] = -dx * hatree_bohr2ev_ang
        f[:, 1] = -dy * hatree_bohr2ev_ang
        f[:, 2] = -dz * hatree_bohr2ev_ang
        e = e * hatree2ev
        self.results = {'energy': e, 'forces': f}


class O4SingletPES(BasePES):
    __pes__func__ = o4_singlet.pot
    """
    O4_singlet
    Y. Paukku, K. R. Yang, Z. Varga, G. Song, J. D. Bender,D. G. Truhlar J.Chem.Phys.
    """
    implemented_properties = [
        "energy",
        "forces", ]


class N4singletPES(BasePES):
    """
System	N4
ID	PES_N4_singlet_umn_v3
Number of bodies	4
Common name	N4(adiabatic ground state)
Functional form	permutation-invariant polynomials
Interface	Section-2
Number of electronic surfaces	1
Number of derivatives	1
References
Bender, Valentini, Nompelis, Paukku, Varga, Truhlar,
Schwartzentruber, and Candler, Journal of Chemical
Physics, submitted. Manuscript no. A15.02.0237
Notes
FORTRAN90 version of n4pes-gpip-meg.f
PES of N4 with special emphasize for
N2 + N2 --> N2 + N + N
- New fit based on extended dataset

    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = n4_singlet.pes_n4_singlet_umn_v3


class O4TripletPES(BasePES):
    """
https://comp.chem.umn.edu/potlib/showPotential.cgi?id=O4_triplet_v2
Potential energy surface information
System	O4
ID	O4_triplet_v2
Number of bodies	4
Common name	O4 triplet (adiabatic ground state)
Functional form	permutation-invariant polynomials
Interface	Section-2
Number of electronic surfaces	1
Number of derivatives	1
References
Y. Paukku, Z. Varga, D. G. Truhlar, to be published (Feb. 09, 2018)
Notes
PES of quintet O4 with special emphasize for
O2 + O2 --> O2 + O + O
- Fit is based on extended dataset
  containing 10180 points
- Mixed-exponential- gaussian (MEG)
  variables are applied to describe the
  long-range interactions
- D3 dispersion correction is implemented for diatoms

    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = o4_triplet.pot


class N2O2tripletPES(BasePES):
    """
    https://comp.chem.umn.edu/potlib/showPotential.cgi?id=N2O2_3A_MB-PIP-MEG2
    O,O,N,N as order of atoms
    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = n2o2_triplet.n2o2_3a_mb_pip_meg2.pot
    __atomic_numbers__ = [8, 8, 7, 7]


if __name__ == '__main__':
    atoms = Atoms('N4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    atoms.calc = N4singletPES()
    atoms.get_potential_energy()
