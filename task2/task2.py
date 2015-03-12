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
from functions import *
from datetime import datetime

time = datetime.now().time()
print(time)

#initialize variables
A = 1000
lmbda = 3
D = 5
mu = 0.5

p = (A, lmbda, D, 2*mu)

print(leastsq(errfunc, p))
print("w for 1 lat 5 eng 0.1")  
time = datetime.now().time()
print(time)
