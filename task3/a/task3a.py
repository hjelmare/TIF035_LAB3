from scipy.optimize import leastsq
from ase import *
from ase.structure import molecule 
from ase.io import read,write
from ase.calculators.eam import EAM
import numpy as np
from math import *
from eam_calculator import get_calc
from datetime import datetime


time = datetime.now().time()
print(time)

#initialize variables
A = 2011
lmbda = 3.21
D = 5.29
mu = 0.93/2

a = 4.05
b = a / 2

p = (A, lmbda, D, 2*mu)
calc = get_calc(p)

#N = 10

for N in range(1,10):
  base = Atoms('Al4', positions=[[0,0,0],[0,b,b],[b,0,b],[b,b,0]])
  mol = Atoms()
  for z in range(N):
    for y in range(N):
      for x in range(N):
        temp = base.copy()
        temp.translate([x*a,y*a,z*a])
        for i in range(4):
          mol.append(temp[i])
  
  mol.set_pbc((True, True, True))
  mol.set_calculator(calc)
  preE = mol.get_potential_energy()

  del mol[0]

  postE = mol.get_potential_energy()
  deltaE = postE - preE
  print(deltaE)

time = datetime.now().time()
print(time)
