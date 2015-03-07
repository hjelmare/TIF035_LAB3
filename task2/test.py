import numpy as np
from scipy.optimize import leastsq



def func(a,x):
  return(a*x**2)

def errfunc(a):
  y = [50.1, 17.9, 32.2]
  x = [5,3,4]
  ret = [0,0,0]
  for i in range(3):
    ret[i] = float(y[i])-float(func(a,x[i]))
  return(ret)


print(leastsq(errfunc, 1 ,))
