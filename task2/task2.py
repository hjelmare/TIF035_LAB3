
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
from ase.calculators.eam import EAM
import numpy as np
from math import *
from gpaw import GPAW
from gpaw import PW

#mol90 = read('POSCAR_0.9')
#  calc = GPAW(mode=PW(50), h= 0.2, xc='PBE', nbands=162, kpts=(int(k),int(k),int(k)), txt='mol90.txt')
#  mol90.set_calculator(calc)

refMol_90 = read('res_POSCAR_0.9.traj');

#initialize variables
A = 10000
lmbda = 2
D = 5
mu = 0.5
rMax = 6
dr=0.2

def Psi(rMax, dr):
  r_c = 6.5
  matrixSize = rMax/dr  

  i = 0
  psi = np.zeros(matrixSize)
  for r in range(0, rMax+dr, dr):
    x = (r-r_c)/3
    if x<=0:
      psi[i] = pow(x,4)/(1+pow(x,4))
      i += 1
    else:
      psi[i] = 0
      i += 1
  return(psi)

def Phi(rMax, dr ,A ,lmbda, psi):
  matrixSize = rMax/dr  
  phi = np.zeros(matrixSize)
  i = 0

  for r in range(0, rMax+dr, dr):
    phi[i] = A*exp(-lmbda*r)
    i += 1
  phi = phi*psi
  return(phi)

def F(D, rho):
  f = -D*np.sqrt(rho)
  return(f)

def Rho(rMax, dr, mu, psi):
  nbrSteps = rMax/dr
  r = np.linspace(0, rMax, nbrSteps)
  rho = np.exp(-2*mu*r)
  rho = rho*psi
  return(rho)

#def optimizationFunc(refMol_90, A,lmbda,D,mu):
#  #test
#  cutoff = 6.28721
#  n = 21
#  rs = np.arrange(0,n)*(cutoff/n)
#  for
psi = Psi(rMax, dr)
phi = Phi(rMax, dr, A, lmbda, psi)
rho = Rho(rMax, dr, mu, psi)
f = F(D, rho)

calc = EAM(elements = mol.get_chemical_symbols(), embedded_energy = f, electron_density = rho, phi=phi)


