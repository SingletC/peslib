import os
from abc import abstractmethod
from typing import Callable, Tuple, List, Optional

import numpy as np
from ase import Atoms
from ase.units import Bohr, Angstrom, Hartree, eV

from peslib.base import BasePES, BasePESv1, DiabaticPES, EvalSurfIO
from peslibf import ch4oh, o4_singlet, n4_singlet, o4_triplet, n2o2_triplet, h2o2, phoh, phsch3, oh3, nh3
from ase.calculators.calculator import Calculator, all_changes
from pathlib import Path

ang2bohr = Angstrom / Bohr
hatree2ev = Hartree / eV
hatree_bohr2ev_ang = hatree2ev * ang2bohr


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


class CH4OH(BasePESv1):
    """
      J. Espinosa-Garcia and J. C. Corchado
      J. Chem. Phys., Vol. 112, p. 5731, 2000
    """
    implemented_properties = [
        "energy",
        "forces", ]
    example_molecule = Atoms('HCH3OH', positions=np.array([[1.59050433, -0.42711963, 2.16479382],
                                                           [0.72147150, 0.12403254, 1.87202512],
                                                           [-0.15513447, -0.44718081, 2.09536767],
                                                           [0.76099165, 0.32025874, 0.82102483],
                                                           [0.68957837, 1.04966117, 2.40691996],
                                                           [0.64375641, 2.29705194, 3.12571355],
                                                           [1.42653445, 2.81342572, 2.92081063], ]

                                                          ))
    __pes__class__ = ch4oh
    __atomic_numbers__ = example_molecule.get_atomic_numbers()


class H2O2(BasePES):
    """
    https://comp.chem.umn.edu/potlib/showPotential.cgi?id=h2o2
    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = h2o2.surf
    __atomic_numbers__ = [1, 8, 8, 1]
    example_molecule = Atoms('HOOH', positions=np.array([[0.839547, 0.880752, 0.422001],
                                                         [0., 0.734058, -0.05275],
                                                         [0., -0.734058, -0.05275],
                                                         [-0.839547, -0.880752, 0.422001]]))

    def _call_method(self, atoms):
        r = atoms.get_positions().flatten() * ang2bohr
        e, f = self.__pes__func__(r)
        f = -f * hatree_bohr2ev_ang
        e = e * hatree2ev
        return e, f.reshape(4, 3)


class PhOH(BasePES):
    """
    only  adiabatic gs for now
    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = phoh.pot
    example_molecule = Atoms('C6OH6',
                             positions=np.array([[-4.58389260, 1.52348991, 0.00000000],
                                                 [-3.18873260, 1.52348991, 0.00000000],
                                                 [-2.49119460, 2.73124091, 0.00000000],
                                                 [-3.18884860, 3.93974991, -0.00119900],
                                                 [-4.58367360, 3.93967191, -0.00167800],
                                                 [-5.28127460, 2.73146591, -0.00068200],
                                                 [-5.29883258, 0.28503907, 0.00058521],
                                                 [-2.63922460, 0.57097691, 0.00131500],
                                                 [-1.39151460, 2.73132091, 0.00063400],
                                                 [-2.63864860, 4.89189291, -0.00125800],
                                                 [-5.13379560, 4.89195291, -0.00263100],
                                                 [-6.38087860, 2.73164891, -0.00086200],
                                                 [-5.45948859, 0.00819239, 0.90565181]]))
    __atomic_numbers__ = example_molecule.get_atomic_numbers()

    def _call_method(self, atoms):
        r = atoms.get_positions() * ang2bohr
        uu, guu, vv, gvv, dvec, cc = self.__pes__func__(1, r.T, 1)
        # for now let us just test gs
        e = vv[0]
        f = - gvv[:, :, 0].T * hatree_bohr2ev_ang
        return e, f


class PhSCH3(BasePES):
    """
    only  adiabatic gs for now
    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = phsch3.pot
    example_molecule = Atoms('C6H5SCH3',
                             positions=np.array([[0.88692142, 1.30720469, 3.58194418],
                                                 [1.52412154, 0.65733868, 2.51631191],
                                                 [2.20802238, 1.40470952, 1.54799459],
                                                 [2.25472475, 2.80194623, 1.64531061],
                                                 [1.61752441, 3.45181226, 2.71094274],
                                                 [0.93362294, 2.70444147, 3.67925965],
                                                 [0.36474787, 0.73657047, 4.32127589],
                                                 [1.48846314, -0.40948255, 2.44200890],
                                                 [2.77689877, 3.37258041, 0.90597921],
                                                 [1.65318199, 4.51863356, 2.78524521],
                                                 [0.44710574, 3.20062865, 4.49289304],
                                                 [3.01736537, 0.57927691, 0.19447144],
                                                 [3.02122618, 1.64338030, -1.23243748],
                                                 [2.01342382, 1.87507482, -1.50730483],
                                                 [3.54483419, 2.54753046, -1.00168405],
                                                 [3.50774135, 1.14719328, -2.04607219]]
                                                ))
    __atomic_numbers__ = example_molecule.get_atomic_numbers()

    def _call_method(self, atoms):
        r = atoms.get_positions() * ang2bohr
        uu, guu, vv, gvv, dvec, cc = self.__pes__func__(1, r.T.astype(np.float64), 1)
        # for now let us just test gs
        e = vv[0]
        f = - gvv[:, :, 0].T * hatree_bohr2ev_ang
        return e, f


class OH3(DiabaticPES):
    """
    OH3 multi state surface
    "Semiclassical Trajectory Studies of Reactive and Nonreactive Scattering of OH(A2Σ+) by H2 Based
    on an Improved Full-Dimensional Ab Initio Diabatic Potential Energy Matrix"
    Shanyu Han, Antonio Gustavo Sampaio de Oliveira Filho, Yinan Shu, Donald G. Truhlar, and Hua Guo
    ChemPhysChem 2022, 23, e202200039
    """
    implemented_properties = [
        "energy",
        "forces", ]
    __pes__func__ = oh3.oh3_pes_truhlar
    example_molecule = Atoms('OH3',
                             positions=np.array([[-1.68972339, -0.33596837, 0.00000000],
                                                 [-0.73375066, -0.42380993, 0.00000000],
                                                 [-2.09163674, 0.53584907, 0.00000000],
                                                 [-2.24378277, -1.11994427, 0.00000000],
                                                 ]
                                                ))

    __atomic_numbers__ = example_molecule.get_atomic_numbers()
    states = 3

    def _call_method(self, atoms):
        r = atoms.get_positions() * ang2bohr
        u, ga = self.__pes__func__(r.T.astype(np.float64))
        # for now let us just test gs

        return u[self.state, self.state], ga[self.state].T * hatree_bohr2ev_ang


class CH2OH(EvalSurfIO):
    """
    CH2OH multi state surface
    "Semiclassical Trajectory Studies of Reactive and Nonreactive Scattering of OH(A2Σ+) by H2 Based
    on an Improved Full-Dimensional Ab Initio Diabatic Potential Energy Matrix"
    Shanyu Han, Antonio Gustavo Sampaio de Oliveira Filho, Yinan Shu, Donald G. Truhlar, and Hua Guo
    ChemPhysChem 2022, 23, e202200039
    """
    implemented_properties = [
        "energy",
        "forces", ]
    example_molecule = Atoms('OCH3',
                             positions=np.array([[1.11580178, 0.49794644, -0.10829221],
                                                 [-1.42643703, 0.05707492, -0.12098057],
                                                 [-2.06853060, -1.86755627, -0.24294965],
                                                 [-2.48213189, 1.43163121, 0.92858405],
                                                 [1.95164196, -0.88179897, -0.94048539], ]
                                                ))

    __atomic_numbers__ = example_molecule.get_atomic_numbers()
    states = 3
    data_dir = Path(__file__).parent.parent / 'src/CH2OH/data/'


class NH3(DiabaticPES):
    """
    OH3 multi state surface
      1) Nangia, S.; Truhlar, D.G. J. Chem. Phys. 2006, 124, 124309
      2) Li, Z.H.; Valero, R., Truhlar, D.G. Theor. Chem. Acc. 2007,
         118, 9-24.
    """
    implemented_properties = [
        "energy",
        "forces", ]
    example_molecule = Atoms('NH3', positions=np.array([[1.59050433, -0.42711963, 2.16479382],
                                                        [0.72147150, 0.12403254, 1.87202512],
                                                        [-0.15513447, -0.44718081, 2.09536767],
                                                        [0.76099165, 0.32025874, 0.82102483]]))
    __pes__func__ = nh3.pot
    __atomic_numbers__ = example_molecule.get_atomic_numbers()
    states = 2

    def _call_method(self, atoms):
        r = atoms.get_positions() * ang2bohr
        u11, u22, u12, v1, v2, gu11, gu22, gu12, gv1, gv2 = self.__pes__func__(r.flatten())
        r = [(v1 * hatree2ev, -gv1.reshape(-1, 3) * hatree_bohr2ev_ang),
             (v2 * hatree2ev, -gv2.reshape(-1, 3) * hatree_bohr2ev_ang)]
        return r[self.state]


if __name__ == '__main__':  # debug purposes
    atoms = CH2OH.example_molecule
    atoms.calc = CH2OH(state=0)
    atoms.get_potential_energy()
    atoms.get_forces()
    from ase.optimize import BFGS

    opt = BFGS(atoms)
    opt.run(fmax=0.01)
