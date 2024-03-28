from unittest import TestCase
from ase.atoms import Atoms

from peslib.pes import O4SingletPES, N4singletPES, O4TripletPES, N2O2tripletPES, CH4OH, H2O2, PhOH, PhSCH3, NH3, OH3


class TestPES(TestCase):
    def test_o4(self):
        atoms = Atoms('O4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
        atoms.calc = O4SingletPES()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_o4_t(self):
        atoms = Atoms('O4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
        atoms.calc = O4TripletPES()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_n4(self):
        atoms = Atoms('N4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
        atoms.calc = N4singletPES()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_n2o2(self):
        atoms = Atoms('O2N2', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.5, 0.5, 0.5)])
        atoms.calc = N2O2tripletPES()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_ch4oh(self):
        atoms = CH4OH.example_molecule
        atoms.calc = CH4OH()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_h2o2(self):
        atoms = H2O2.example_molecule
        atoms.calc = H2O2()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_phoh(self):
        atoms = PhOH.example_molecule
        atoms.calc = PhOH()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')

    def test_oh3(self):
        for i in [0, 1, 2]:  # three states
            atoms = OH3.example_molecule
            atoms.calc = OH3(state=i)
            print(f'E : {atoms.get_potential_energy()}')
            print(f'Force : {atoms.get_forces()}')

    def test_nh3(self):
        for i in [0, 1]:  # two states
            atoms = NH3.example_molecule
            atoms.calc = NH3(state=i)
            print(f'E : {atoms.get_potential_energy()}')
            print(f'Force : {atoms.get_forces()}')

    def test_phsch3(self):
        ...  # somehow this is not working. try gcc 4.8 ?
        # atoms = PhSCH3.example_molecule
        # atoms.calc = PhSCH3()
        # print(f'E : {atoms.get_potential_energy()}')
        # print(f'Force : {atoms.get_forces()}')
