import io
import subprocess
from typing import Callable, Tuple, List, Optional

import numpy as np
from ase import Atoms

from peslib.utils import atoms_to_eval_surf_io, ang2bohr, hatree_bohr2ev_ang, hatree2ev
from ase.calculators.calculator import Calculator, all_changes
from pathlib import Path




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

    def _call_method(self, atoms):
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
        return e, f

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            self.atoms = atoms.copy()
        self.check_atomic_numbers(atoms)
        e, f = self._call_method(atoms)
        self.results = {'energy': e, 'forces': f}


class AdiabaticPES(BasePES):
    __pes__func__: Callable[
        [np.ndarray, np.ndarray, np.ndarray], Tuple[float, np.ndarray, np.ndarray, np.ndarray]] = None
    __atomic_numbers__: Optional[List[int]] = None
    example_molecule: Optional[Atoms] = None
    states: int = 0

    def __init__(self, state: int = 0, **kwargs):
        Calculator.__init__(self, **kwargs)
        if state >= self.states:
            raise ValueError(f'{self.__class__.__name__} only support {self.states} states')
        self.state = state

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        if atoms is not None:
            self.atoms = atoms.copy()
        self.check_atomic_numbers(atoms)
        e, f, u, u_grad = self._call_method(atoms)
        self.results = {'energy': e, 'forces': f[self.state,self.state],'nacdr':f,
                        'u':u, 'u_grad':u_grad}


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


class EvalSurfIO(AdiabaticPES):
    """
    EvalSurf Type PES call external binary file read data and output multiple states PES
    """

    data_dir: Path
    executable_path = Path(__file__).parent / 'evalsurf.x'

    def _call_method(self, atoms):
        # create in-memory tmp file
        in_memory_file = io.StringIO()
        input = atoms_to_eval_surf_io(atoms)
        result = subprocess.run([self.executable_path], input=input, text=True, capture_output=True,
                                cwd=self.data_dir)
        if result.stderr:
            raise PESLIBError(f'Error in EvalSurf: {result.stderr}')

        output = result.stdout
        # parse output
        lines = output.split('\n')
        v_lines = lines.index(' Adiabatic energy(a.u.)')
        if v_lines == -1:
            raise PESLIBError('Adiabatic energy not found')
        v = np.fromstring(lines[v_lines + 1], dtype=float, sep=' ')
        u_lines = lines.index(' Quasi-diabatic Hamiltonian')
        if u_lines == -1:
            raise PESLIBError('Quasi-diabatic Hamiltonian not found')
        u = np.array([[float(val) for val in line.split()] for line in lines[u_lines + 1:u_lines+self.states+1]])
        grad = np.empty((len(v),len(v), len(atoms), 3))
        u_grad = np.empty((len(u),len(u), len(atoms), 3))
        for i in lines[u_lines:]:
            if 'Adiabatic gradients of states:' in i:
                gau_ls = []
                n_state, m_state = np.fromstring(i[-5:], dtype=int, sep=' ')
                i_idx = lines.index(i)
                while ' ' != lines[i_idx+1]:
                    i_idx += 1
                    gau_ls.append(np.fromstring(lines[i_idx], dtype=float, sep=' '))
                if not gau_ls:
                    raise PESLIBError('Adiabatic gradients not found')
                gv = np.array(gau_ls)
                if n_state!=m_state:
                    gv = - gv * ang2bohr # to Ang**-1
                else:
                    gv = gv  * hatree_bohr2ev_ang
                grad[n_state-1, m_state-1] = gv
                grad[m_state-1, n_state-1] = gv
            if 'Diabatic gradients of states' in i:
                gu_ls = []
                n_state, m_state = np.fromstring(i[-5:], dtype=int, sep=' ')
                i_idx = lines.index(i)
                while ' ' != lines[i_idx+1]:
                    i_idx += 1
                    gu_ls.append(np.fromstring(lines[i_idx], dtype=float, sep=' '))
                if not gau_ls:
                    raise PESLIBError('diabatic gradients not found')
                gu = np.array(gu_ls)
                gu = gu  * hatree_bohr2ev_ang
                u_grad[n_state-1, m_state-1] = gu
                u_grad[m_state-1, n_state-1] = gu

        return v[self.state]*hatree2ev,  - grad, u,u_grad # force is negative gradient


class PESLIBError(Exception):
    pass