from scipy.optimize import leastsq
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
from ase.calculators.eam import EAM
import numpy as np
from math import *
from gpaw import GPAW
from gpaw import PW
from eam_calculator import get_calc

# we need to shape the force vectors to the appropriate form
# and maybe think about other properties as well

def getReferenceValues():
  mol = read('res_POSCAR_1.0.traj')
  forces = mol.get_forces()
  forces = forces.reshape(1,forces.size)
  # later, we need to load the other files as well...  
  return(forces)

def getEamValues(p):
  mol = read('res_POSCAR_1.0.traj')
  calc = get_calc(p)
  mol.set_calculator(calc)
  forces = mol.get_forces()
  forces = forces.reshape(1,forces.size)
  
  return(forces)


def errfunc(p):
  Aref = getReferenceValues()
  Aref = Aref[0]
  Aeam = getEamValues(p)
  Aeam = Aeam[0]
  ret = [0]
  for i in range(Aref.size):
    ret = ret + [0]

  for i in range(Aref.size):
    ret[i] = float(Aref[i])-float(Aeam[i])

  return( ret )

# end of errfunc declaration

#initialize variables
A = 10000
lmbda = 2
D = 5
mu = 0.5

p = (A, lmbda, D, 2*mu)

print(leastsq(errfunc, p))
  
