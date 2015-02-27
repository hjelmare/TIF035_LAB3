from scipy.optimize import leastsq



def func(a,x):
  return(a[0]*x + a[1])

def errfunc(a,x,y):
  return(y-func(a,x))

print(leastsq(errfunc, [1,2] , args=([1, 2, 3],[ 1, 3, 5])))
