from ase import *
#from ase.structure import molecule 
from ase.io import read,write
from ase.calculators.eam import EAM
import numpy as np
from math import *
from eam_calculator import get_calc
from datetime import datetime
from ase.neb import NEB
from ase.optimize import MDMin
from ase.optimize import QuasiNewton

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
#calc = GPAW(mode=PW(200))

N = 10

#initial = read('res_POSCAR_1.0.traj')
start = Atoms('Al4', positions=[[0,0,0],[0,b,b],[b,0,b],[b,b,0]])
initial = Atoms()
for z in range(N):
  for y in range(N):
    for x in range(N):
      temp = start.copy()
      temp.translate([x*a,y*a,z*a])
      for i in range(4):
        print((x,y,z,i))
        initial.append(temp[i])

tempPos = initial[0].position
del initial[0]
final = initial.copy()
final[1].position = tempPos

#QuasiNewton(initial).run(fmax=0.05)
#QuasiNewton(final).run(fmax=0.05)

images = [initial]
images += [initial.copy() for i in range(3)]
images += [final]

neb = NEB(images)

neb.interpolate()

for image in images[1:4]:
  image.set_calculator(calc)

print("its go time")
optimizer = MDMin(neb, trajectory='neb.traj')
optimizer.run(fmax=0.04)

time = datetime.now().time()
print(time)
