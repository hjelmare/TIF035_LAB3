
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule #kanske inte behövs
from ase.io import read,write
import nmpy as np
from gpaw import GPAW
from gpaw import PW

mol90 = read('POSCAR_0.9')

print(mol90.positions)




