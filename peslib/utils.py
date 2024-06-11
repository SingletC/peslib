from ase import Atoms, Atom
from ase.units import Bohr, Angstrom, Hartree, eV
import numpy as np
from ase import Atoms
try:
    # risky import for development
    from scipy.optimize._numdiff import approx_derivative
except ImportError:
    pass

ang2bohr = Angstrom / Bohr
hatree2ev = Hartree / eV
hatree_bohr2ev_ang = hatree2ev * ang2bohr


def atoms_to_eval_surf_io(atoms: Atoms) -> str:
    """
    Convert ASE Atoms to EvalSurf input file
    """
    r = []
    atoms_ = atoms.copy()
    atoms_.positions *= ang2bohr
    for i in atoms_:
        i: Atom
        r.append(f' {i.symbol}      {i.number:1.0f}.{i.position[0]:14.8f}{i.position[1]:14.8f}{i.position[2]:14.8f}'
                 f'{i.mass:14.8f}')
    return '\n'.join(r)


def num_gradient(atoms: Atoms):
    """
    calculate numerical gradient of Ase Atoms with calculator
    """
    old_pos = atoms.positions.copy()
    def get_e(pos: np.array):
        atoms.positions = pos.reshape(-1, 3)
        return atoms.get_potential_energy()

    grad = approx_derivative(get_e, atoms.positions.flatten(), method='3-point', abs_step=1e-5).reshape(-1, 3)
    atoms.positions = old_pos
    return grad
