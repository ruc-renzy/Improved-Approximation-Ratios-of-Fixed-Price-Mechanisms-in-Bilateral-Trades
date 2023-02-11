import numpy as np
import math
import time


beta = 0.7
g1 = (1-beta)/beta
n_e = 50
n_t = 50
sol = np.zeros((n_e*n_t+1, n_e*n_t+1, n_e+1))

x = np.linspace(0, 1, n_e*n_t+1)
y = np.linspace(0, 1, n_e*n_t+1)
z = np.linspace(0, 1, n_e+1)

z_ = np.linspace(1/n_t, 1, n_t)


def cal_fp(hs, ks, gs):
    Hs = hs/n_e
    Ks = ks/n_e/n_t
    Gs = gs/n_e/n_t + g1
    numer = (1-beta)*g1-beta*Ks+beta*Hs*Gs
    return numer/(Gs*(Gs+Hs/n_t))/n_t

t_s = time.time()
for i in z_:
    sol_t = np.zeros((n_e*n_t+1, n_e*n_t+1, n_e+1))
    for x_index in range(n_e*n_t+1):
        for y_index in range(n_e*n_t+1):
            for z_index in range(n_e+1):
                if (n_e*n_t*x_index+n_e+n_t)*(1-i)<y_index**2 or (y_index+n_e+n_t)*z_index<n_e*x_index:
                    continue
                else:
                    for hs in range(z_index+1):
                        ks = x_index - int(hs**2/n_e)
                        gs = y_index - hs
                        if ks >= 0 and gs >= 0:
                            temp = cal_fp(hs, ks, gs)
                            if temp + sol[x_index][y_index][z_index] > sol_t[ks][gs][hs]:
                                sol_t[ks][gs][hs] = temp + sol[x_index][y_index][z_index]
    sol = sol_t

max_fp = 0.0
for x_index in range(n_e*n_t + 1):
    for y_index in range(n_e*n_t + 1):
        for z_index in range(n_e + 1):
            if sol[x_index][y_index][z_index] > max_fp:
                max_fp = sol[x_index][y_index][z_index]

print(time.time()-t_s)
print(beta)
print(n_e, n_t)
print(max_fp)
