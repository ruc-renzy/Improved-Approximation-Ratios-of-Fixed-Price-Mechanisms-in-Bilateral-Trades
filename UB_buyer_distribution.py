import math
from scipy import integrate


r = 6
b0 = 0.25
b1 = 0.6
lam = 1/(1-math.e**(r*(b0-b1)))

def cal_HB(x):
    return 1 - lam*(1-math.e**(r*(b0-x)))

def cal_GB(x):
    F = lambda y: 1 - lam * (1 - math.e ** (r * (b0 - y)))
    v, err = integrate.quad(F, x, b1)
    return v

def cal_HB2(x):
    F = lambda y: (1 - lam*(1-math.e**(r*(b0-y))))**2
    v, err = integrate.quad(F, x, b1)
    return v


beta = 0.7381
G1 = (1-beta)/beta
v1 = (1-b1)*beta
print('first part:', v1)


F2 = lambda x:beta*(cal_HB(x)/(cal_GB(x)+G1)-cal_HB2(x)/(cal_GB(x)+G1)**2)+(1-beta)*G1/(cal_GB(x)+G1)**2
v2, err2 = integrate.quad(F2, b0, b1)

print('second part:', v2)


A = beta*(cal_GB(b0)+G1-cal_HB2(b0))+(1-beta)*G1
v3 = b0*A/(b0*(cal_GB(b0)+G1)+(cal_GB(b0)+G1)**2)
print('third part:', v3)

print('total integral:', v1+v2+v3)

