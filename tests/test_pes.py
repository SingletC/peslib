from unittest import TestCase
from ase.atoms import Atoms

from peslib.pes import O4SingletPES


class TestO4PES(TestCase):
    def test_calculate(self):
        atoms = Atoms('O4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
        atoms.calc = O4SingletPES()
        print(f'E : {atoms.get_potential_energy()}')
        print(f'Force : {atoms.get_forces()}')
