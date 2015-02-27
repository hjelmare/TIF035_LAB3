
from ase import *
from ase.optimize import BFGS
from ase.structure import molecule 
from ase.io import read,write
from ase.calculators.eam import EAM
import numpy as np
from math import *
from gpaw import GPAW
from gpaw import PW
from func_task2 import *

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

# here's where we had the function declarations earlier

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


