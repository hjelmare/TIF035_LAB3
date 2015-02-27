
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
import numpy as np
from gpaw import GPAW
from gpaw import PW

mol90 = read('POSCAR_0.9')

for k in (1, 2, 3):
  calc = GPAW(mode=PW(50), h= 0.2, xc='PBE', nbands=162, kpts=(int(k),int(k),int(k)), txt='mol90.txt')

  mol90.set_calculator(calc)
  print(mol90.get_forces())



