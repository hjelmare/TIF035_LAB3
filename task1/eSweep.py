
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
import numpy as np
from gpaw import GPAW
from gpaw import PW

mol90 = read('POSCAR_0.9')

energy = np.linspace(50, 200, 7)

for e in energy:
  calc = GPAW(mode=PW(int(e)), h= gridWidth, xc='PBE', nbands=162, kpts=(1,1,1), txt='mol90_eSweep.txt')

  mol90.set_calculator(calc)
  print(mol90.get_forces())



