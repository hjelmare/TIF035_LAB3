
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
  # later, we need to load the other files as well...  
  return(forces)

def getEamValues(p):
  mol = read('res_POSCAR_1.0.traj')
  calc = get_calc(p)
  mol.set_calculator(calc)
  forces = mol.get_forces()
  
  return(forces)


def errfunc(p):
  Aref = getReferenceValues()

  Aeam = getEamValues()  

  return( Aref-Aeam )

# end of errfunc declaration

#mol90 = read('POSCAR_0.9')
#  calc = GPAW(mode=PW(50), h= 0.2, xc='PBE', nbands=162, kpts=(int(k),int(k),int(k)), txt='mol90.txt')
#  mol90.set_calculator(calc)

refMol_90 = read('res_POSCAR_0.9.traj');

#initialize variables
A = 10000
lmbda = 2
D = 5
mu = 0.5

p = (A, lmbda, D, 2*mu)
calc = get_calc(p)

mol.set_calculator(calc)

  
