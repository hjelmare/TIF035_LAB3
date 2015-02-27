
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
import numpy as np
from math import *
from gpaw import GPAW
from gpaw import PW

#mol90 = read('POSCAR_0.9')
#  calc = GPAW(mode=PW(50), h= 0.2, xc='PBE', nbands=162, kpts=(int(k),int(k),int(k)), txt='mol90.txt')
#  mol90.set_calculator(calc)

refMol_90 = read('res_POSCAR_0.9.traj');

def psi(r):
  r_c = 6.5
  x = (r-r_c)/3
  if x<=0
    return(math.pow(x,4)/(1+math.pow(x,4)))
  else
    return(0)

def Phi(r,A,lmbda,psi):
  return(A*exp(-lmbda*r)*psi)

def F(D, rho):
  return(-D*math.sqrt(rho))

def Rho(mu, r, psi):
  return(exp(-2*mu*r)*psi)

print(psi(2))

#def optimizationFunc(refMol_90, A,lmbda,D,mu)





