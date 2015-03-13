from scipy.optimize import leastsq
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
from ase.calculators.eam import EAM
from ase.lattice import bulk
import numpy as np
from math import *
from gpaw import GPAW
from gpaw import PW
from eam_calculator import get_calc
from functions import *
from datetime import datetime

time = datetime.now().time()
print(time)

latticeParam = 4.05

#initialize variables
A = 2011
lmbda = 3.21
D = 5.29
mu = 0.93/2

p = (A, lmbda, D, 2*mu)

dE1m2 = []
dE4 = []
dE1p2 = []
e = []
for deformation in (0.2,0.1,0.001):
  mol = read('res_POSCAR_1.0.traj')
#  mol = bulk('Al','fcc',a = latticeParam)

  calc = get_calc(p)
  mol.set_calculator(calc)

  preEnergy = mol.get_potential_energy()
  preVolume = mol.get_volume()

  strain_change1m2 = np.array(((deformation, 0, 0),(0,-deformation, 0),(0,0,deformation**2/(1-deformation**2))))
  strain_change4 = np.array(((0,0.5*deformation, 0),(0.5*deformation,0, 0),(0,0,deformation**2/(4-deformation**2))))
  strain_change1p2 = np.array(((deformation, 0, 0),(0,deformation, 0),(0,0,0)))

  strain1m2 = np.identity(3) + strain_change1m2
  strain4 = np.identity(3) + strain_change4
  strain1p2 = np.identity(3) + strain_change1p2
  
  preCell = mol.get_cell()
  
  mol.set_cell(np.dot(strain4, preCell), scale_atoms=True)
  postEnergy4 = mol.get_potential_energy()
  
  mol.set_cell(np.dot(strain1m2, preCell), scale_atoms=True)
  postEnergy1m2 = mol.get_potential_energy()
  
  mol.set_cell(np.dot(strain1p2, preCell), scale_atoms=True)
  postEnergy1p2 = mol.get_potential_energy()

  dE1m2 += [postEnergy1m2 - preEnergy]
  dE4 += [postEnergy4 - preEnergy]
  dE1p2 += [postEnergy1p2 - preEnergy]
  e += [deformation]

C1m2 = np.polyfit(e,dE1m2,2)[0]
C4 =  np.polyfit(e,dE4,2)[0]
C1p2 = np.polyfit(e,dE1p2,2)[0]

conversion = 1.66e11

c12 = conversion*(C1p2 - C1m2)/2
c11 = conversion*C1m2+c12
c44 = conversion*C4


print((c11, c12, c44))

B = c11+2*c12 / 3
tetra = (c11-c12)/2

print((B, tetra))


time = datetime.now().time()
print(time)
