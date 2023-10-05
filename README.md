## ASE calculator type PES surface
repo currently under development
current supported PES:
- [O4 singlet](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=O4_singlet)
- [O4 triplet](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=O4_triplet_v2)
- [N4_singlet](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=PES_N4_singlet_umn_v3)
- [N2O2_triplet](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=PES_N2O2_triplet_umn_v3)
- [CH4OH](https://doi.org/10.1063/1.481148)
- [H2O2](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=h2o2)
- [PhOH](https://comp.chem.umn.edu/potlib/showPotential.cgi?id=phoh_aprp) *adiabatic gs
## Installation
make sure there is gfortran compiler installed
```bash
gfortran --version
```
if not, install it, then
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