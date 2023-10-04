from typing import Callable, Tuple, List, Optional

import numpy as np
from ase import Atoms
from ase.units import Bohr, Angstrom, Hartree, eV
from peslibf import ch4oh, o4_singlet, n4_singlet, o4_triplet, n2o2_triplet,c2h6oh
from ase.calculators.calculator import Calculator, all_changes

ang2bohr = Angstrom / Bohr
hatree2ev = Hartree / eV
hatree_bohr2ev_ang = hatree2ev * ang2bohr


class BasePES(Calculator):
    __pes__func__: Callable[
        [np.ndarray, np.ndarray, np.ndarray], Tuple[float, np.ndarray, np.ndarray, np.ndarray]] = None
    __atomic_numbers__: Optional[List[int]] = None
    example_molecule: Optional[Atoms] = None

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def check_atomic_numbers(self, atoms: Atoms):
        if self.__atomic_numbers__ is None:
            return
        if not (atoms.get_atomic_numbers() == self.__atomic_numbers__).all():
            raise ValueError(
                f'order of Atomic numbers of atoms should be {self.__atomic_numbers__}')

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            self.atoms = atoms.copy()
        self.check_atomic_numbers(atoms)
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


class BasePESv1(BasePES):
    """
    POTLIB type 1 interface.
    """
    __pes__class__: object = None
    _NASURF = np.zeros((6, 6), dtype=int)
    _NASURF[0] = 1
    _NFLAG = np.zeros(20).astype(int)
    _NFLAG[0] = 1
    _NFLAG[1] = 1
    _NDER = 1  # 1st gradient

    def __init__(self):
        BasePES.__init__(self)
        self.__pes__class__.prepot()

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            self.atoms = atoms.copy()
        self.check_atomic_numbers(atoms)
        r = atoms.get_positions().flatten() * ang2bohr
        # make it 75 array by filling zero
        r = np.hstack((r, np.zeros(75 - len(r)))).reshape(-1, 3)

        self.__pes__class__.usricm.cart = r
        self.__pes__class__.usricm.NFLAG = self._NFLAG
        self.__pes__class__.usricm.nder = self._NDER
        self.__pes__class__.usricm.nasurf = self._NASURF

        self.__pes__class__.pot()
        e, dx = self.__pes__class__.usrocm.pengygs, self.__pes__class__.usrocm.dgscart  # why dedx is in the middle?
        dx = dx[~(dx == 0).all(axis=1)]  # dangerous, if some row indeed is zero will be removed.
        f = -dx * hatree_bohr2ev_ang
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


# class C2H6OtripletPES(BasePESv1):
#     """
#     https://comp.chem.umn.edu/potlib/showPotential.cgi?id=c2h6oh-v2020
# C        This potential is written such that:
# C                       X(1)  - X(3)  : X, Y, Z for C1
# C                       X(4)  - X(6)  : X, Y, Z for C2
# C                       X(7)  - X(9)  : X, Y, Z for H1
# C                       X(10) - X(12) : X, Y, Z for H2
# C                       X(13) - X(15) : X, Y, Z for H3
# C                       X(16) - X(18) : X, Y, Z for H4
# C                       X(19) - X(21) : X, Y, Z for H5
# C                       X(22) - X(24) : X, Y, Z for H6
# C                       X(25) - X(27) : X, Y, Z for O
# C                       X(28) - X(30) : X, Y, Z for HO
#     """
#     implemented_properties = [
#         "energy",
#         "forces", ]
#     __pes__func__ = c2h6oh.pot
#     __atomic_numbers__ = [6, 6] + [1] * 6 + [8, 1]
#     c2h6 = molecule('C2H6')
#     x = c2h6.get_positions()
#     x = np.vstack((x, np.array([[0, 0, 2], [0 ,0, 3]])))
#     example_molecule = Atoms('C2H6OH', positions=x)

class CH4OH(BasePESv1):
    """
      J. Espinosa-Garcia and J. C. Corchado
      J. Chem. Phys., Vol. 112, p. 5731, 2000
    """
    implemented_properties = [
        "energy",
        "forces", ]
    example_molecule = Atoms('HCH3OH', positions=np.array([[0., 0., 0.],
                                                           [1.88972613, 0., 0.],
                                                           [0., 1.88972613, 0.],
                                                           [0.94486306, 0.94486306, 0.94486306],
                                                           [0.94486306, 0.94486306, 2.83458919],
                                                           [0.94486306, 0.94486306, 4.72431531],
                                                           [0.94486306, 0.94486306, 6.61404144], ]

                                                          ))
    __pes__class__ = ch4oh
    __atomic_numbers__ = example_molecule.get_atomic_numbers()


class C2H6OH(BasePESv1):
    """
    C. Rangel, M. Garcia-Chamorro, J.C. Corchado and J. Espinosa-Garcia, PCCP, 22 14796 (2020)
      https://comp.chem.umn.edu/potlib/showPotential.cgi?id=c2h6oh-v2020
    """
    implemented_properties = [
        "energy",
        "forces", ]
    example_molecule = Atoms('C2H6OH', positions=np.array([[-2.70469798, 0.12751678, 0.00000000],
                                                           [-2.19135576, 0.85347305, 1.25740497],
                                                           [-2.34804355, -0.88129323, 0.00000000],
                                                           [-2.34802514, 0.63191497, -0.87365150],
                                                           [-3.77469798, 0.12752996, 0.00000000],
                                                           [-1.12135757, 0.85176598, 1.25838332],
                                                           [-2.54641232, 1.86284607, 1.25642758],
                                                           [-2.54962472, 0.35020505, 2.13105517],
                                                           [-2.02369196, 2.85589893, 2.47385774],
                                                           [-1.75811093, 3.34840994, 3.25392159]]

                                                          ))
    __pes__class__ = c2h6oh
    __atomic_numbers__ = example_molecule.get_atomic_numbers()


if __name__ == '__main__':
    atoms = Atoms('N4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
    atoms.calc = N4singletPES()
    atoms.get_potential_energy()
