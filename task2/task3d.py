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
from ase.neb import NEB
from ase.optimize import MDMin
from ase.optimize import QuasiNewton

time = datetime.now().time()
print(time)
print('jobbigt')

#initialize variables
A = 2011
lmbda = 3.21
D = 5.29
mu = 0.93/2

a = 4.05
b = a / 2

p = (A, lmbda, D, 2*mu)
calc = get_calc(p)
calc = GPAW(mode=PW(200))

initial = read('res_POSCAR_1.0.traj')
#initial = Atoms('Al', positions=[[0,0,0],[0,b,b],[b,0,b],[b,b,0]]
tempPos = initial[0].position
del initial[0]
final = initial.copy()
final[1].position = tempPos

QuasiNewton(initial).run(fmax=0.05)
QuasiNewton(final).run(fmax=0.05)

images = [initial]
images += [initial.copy() for i in range(3)]
images += [final]

neb = NEB(images)

neb.interpolate()

for image in images[1:4]:
  image.set_calculator(calc)

#optimizer = BFGS(neb, trajectory='neb.traj')
optimizer = MDMin(neb, trajectory='neb.traj')
optimizer.run(fmax=0.04)

time = datetime.now().time()
print(time)
