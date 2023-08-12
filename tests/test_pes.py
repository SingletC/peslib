from unittest import TestCase
from ase.atoms import Atoms

from peslib.pes import O4SingletPES,N4singletPES,O4TripletPES


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

