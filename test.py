from ase import *
import numpy as np
from ase.io import PickleTrajectory
from ase.io import read, write

mol = read('mol_task1.traj')

print(mol.get_forces())
