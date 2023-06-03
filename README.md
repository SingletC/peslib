## ASE calculator type PES surface
repo currently under development
current supported PES:
- [O4 singlet](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=O4_singlet)
## Installation
make sure there is gfortran compiler installed
```bash
gfortran --version
```
if not, install it
```bash
pip install .
```
verify installation by running
```bash
cd tests  
python -m unittest
```

## Usage

```python
from ase.atoms import Atoms
from peslib.pes import O4SingletPES

atoms = Atoms('O4', positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
atoms.calc = O4SingletPES()
print(f'E : {atoms.get_potential_energy()}')
print(f'Force : {atoms.get_forces()}')
```