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
  
  lattice = mol.get_cell()[0,0] / 3  #for EAM we get the minimizing value, so fix that here as well..........

  energy =  mol.get_potential_energy()
  
  ret = np.concatenate(([[lattice]],[[energy]], forces),1)
  # later, we need to load the other files as well...
  return(ret)



def getEamValues(p):
  mol = read('res_POSCAR_1.0.traj')
  calc = get_calc(p)
  mol.set_calculator(calc)
  forces = mol.get_forces()
  energy = [mol.get_potential_energy()]

  eps = 0.02
  energies = []
  lattice = []
  initialLattice = mol.get_cell()[0,0]
  for scale in np.linspace(1-eps, 1+eps, 3):
    mol.set_cell(3*[initialLattice*scale], scale_atoms=True)
    energies = energies + [ mol.get_potential_energy()]
    lattice = lattice + [[initialLattice*scale/3.0]]
  
  lattice = np.array(lattice)
  energies = np.array(energies)

  functions = np.array([lattice**0, lattice, lattice**2])
  polyfit = np.linalg.lstsq(functions.T[0], energies)[0]
  derivative = [polyfit[1], 2*polyfit[2]]
  a0 = np.linalg.solve([derivative[1]], [-derivative[0]] )
  forces = forces.reshape(1,forces.size)
  print(a0)
  print(energy)

  ret = np.concatenate((a0, energy, forces),1)
  
  return(ret)


def errfunc(p):
  Aref = getReferenceValues()
  Aref = Aref[0]
  Aeam = getEamValues(p)
  Aeam = Aeam[0]
  ret = [0]
  w = [1]
  for i in range(Aref.size):
    ret = ret + [0]
    w = w + [w[0]]
  
  w[1] = 50   # lattice param, remember to change task2py output as well
  w[2] = 1   # energy

  for i in range(Aref.size):
    ret[i] = float(Aref[i])-float(Aeam[i])
    ret[i] = ret[i]*w[i]
  return( ret )
  
