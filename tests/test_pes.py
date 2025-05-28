from unittest import TestCase
from ase.atoms import Atoms
import numpy as np

from peslib.pes import O4SingletPES, N4singletPES, O4TripletPES, N2O2tripletPES, CH4OH, H2O2, PhOH, PhSCH3, NH3, OH3, \
    CH2OH, O2H3


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

    def test_ch2oh(self):
        for i in [0, 1, 2]:  # two states
            atoms = CH2OH.example_molecule
            atoms.calc = CH2OH(state=i)
            print(f'E : {atoms.get_potential_energy()}')
            print(f'Force : {atoms.get_forces()}')

    def test_phsch3(self):
        ...  # somehow this is not working. try gcc 4.8 ?
        # atoms = PhSCH3.example_molecule
        # atoms.calc = PhSCH3()
        # print(f'E : {atoms.get_potential_energy()}')
        # print(f'Force : {atoms.get_forces()}')

    def test_o2h3(self):
        """Test O2H3 PIP-NN PES for OH + H2O reaction system."""
        # Test with example molecule
        atoms = O2H3.example_molecule
        atoms.calc = O2H3()
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        
        print(f'O2H3 Example - E : {energy:.6f} eV')
        print(f'O2H3 Example - Force shape : {forces.shape}')
        print(f'O2H3 Example - Max force magnitude: {np.max(np.abs(forces)):.6f} eV/Å')
        
        # Test with a different geometry (more realistic OH + H2O)
        atoms2 = Atoms('HHHOO', positions=np.array([
            [0.00000000, 0.00000000, 0.00000000],  # H (from OH)
            [1.13322764, 0.00000000, 1.22881848],  # H (from H2O)
            [1.88236977, 1.46895492, 1.50296719],  # H (from H2O)  
            [0.00000000, 0.00000000, 0.96684775],  # O (from OH)
            [2.10278967, 0.62594543, 1.08398697]   # O (from H2O)
        ]))
        
        atoms2.calc = O2H3()
        energy2 = atoms2.get_potential_energy()
        forces2 = atoms2.get_forces()
        
        print(f'O2H3 Realistic - E : {energy2:.6f} eV')
        print(f'O2H3 Realistic - Force shape : {forces2.shape}')
        print(f'O2H3 Realistic - Max force magnitude: {np.max(np.abs(forces2)):.6f} eV/Å')
        
        # Basic validation
        assert isinstance(energy, float)
        assert isinstance(energy2, float)
        assert forces.shape == (5, 3)
        assert forces2.shape == (5, 3)
        assert not np.isnan(energy)
        assert not np.isnan(energy2)
        assert not np.any(np.isnan(forces))
        assert not np.any(np.isnan(forces2))
