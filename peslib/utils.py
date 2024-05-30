from ase import Atoms, Atom
from ase.units import Bohr, Angstrom, Hartree, eV
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
