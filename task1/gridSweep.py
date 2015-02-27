
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
import numpy as np
from gpaw import GPAW
from gpaw import PW

mol90 = read('POSCAR_0.9')


for gridWidth in (0.5, 0.4, 0.3, 0.2, 0.1):
  calc = GPAW(mode=PW(50), h= gridWidth, xc='PBE', nbands=162, kpts=(1,1,1), txt='mol90_gridSweep.txt')

  mol90.set_calculator(calc)
  print(mol90.get_forces())



