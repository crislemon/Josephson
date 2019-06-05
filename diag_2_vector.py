#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 12:13:28 2019

@author: cristina
"""

import numpy as np
from numpy import linalg as LA
diag = LA.eigh
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
import time

pi = np.pi
t1 = time.time()

N = 3000 #number of sites
m = 1.0 #effective mass
delta =1.35/27211.6 #SC gap
mu = 1.0/27211.6 #chemical potential
mu = 0.0
a = 4.98/0.529 ##lattice constant

H = np.zeros([N * 2,N  * 2], dtype=complex)
h = np.zeros([N, N, 2, 2], dtype=complex)

factor = 1/(m*a**2) - mu
factor_2 = -1/(2*m*a**2)

t = np.arange(N)
T = np.meshgrid(t)
t_i = np.reshape(T, (((N)), ))

h[t_i, t_i, 0, 0] = factor
h[t_i, t_i, 1, 1] = - factor

h[t_i, t_i, 0, 1] = delta
h[t_i, t_i, 1, 0] = delta
   

r = np.arange(N-1)
r2 = r + 1
R = np.meshgrid(r)
R2 = np.meshgrid(r2)
r_i = np.reshape(R, (N-1), )
r_j = np.reshape(R2, (N-1), )


h[r_i, r_j, 0, 0] = factor_2 
h[r_i, r_j, 1, 1] = - factor_2 

h[r_j, r_i, 0, 0] = factor_2 
h[r_j, r_i, 1, 1] = - factor_2 


for i in range(N):
        for j in range(N):
            for t_i in range(2):
                for t_j in range(2):
                    H[(i) * 2 + t_i, (j) * 2 + t_j] = h[i, j, t_i, t_j]
                    
H = np.matrix(H)  
T = np.allclose(H, H.getH())###check if Hermitian
print('Is H an Hermitian matrix?', T)

(E, psi) = diag(H)####diagonalize H

u_1 = np.zeros(int(len(E)))
v_1 = np.zeros(int(len(E)))

u_M = np.zeros(int(len(E)))
v_M = np.zeros(int(len(E)))

M = int(N/2.0)

for i in range(int(len(E))):
    
    u_1[i] = psi[0,i]
    v_1[i] = psi[1,i]
    
    u_M[i] = psi[2*(M-1),i]
    v_M[i] = psi[2*M-1,i]


####LDOS functions
def LDOS_up(omega, E, u, Damping):
    t = sum ( u**2 / (omega - E + 1j*Damping) )
    tt = -1/pi*np.imag(t)
    return(tt)
    
def LDOS_down(omega, E, v, Damping):
    t = sum ( v**2 / (omega + E + 1j*Damping) )
    tt = -1/pi*np.imag(t)
    return(tt)
    

###LDOS calculation
omega = np.linspace(-3*delta, 3*delta, 2000)    

LDOS_1_up = np.zeros(len(omega))
LDOS_1_down = np.zeros(len(omega))

LDOS_M_up = np.zeros(len(omega))
LDOS_M_down = np.zeros(len(omega))

D = 0.02/27211.6
for i in range(len(omega)):
    LDOS_1_up[i] = LDOS_up(omega[i], E, u_1, D)    
    LDOS_1_down[i] = LDOS_up(omega[i], E, v_1, D)

    LDOS_M_up[i] = LDOS_up(omega[i], E, u_M, D)    
    LDOS_M_down[i] = LDOS_up(omega[i], E, v_M, D)
    
plt.figure(1)
plt.plot(omega*27211.6, LDOS_1_up + LDOS_1_down) 
plt.plot(omega*27211.6, LDOS_1_up, label = 'up')  
plt.plot(omega*27211.6, LDOS_1_down, label = 'down')
plt.title('Site 1')     
plt.legend()    
    
plt.figure(2)
plt.plot(omega*27211.6, LDOS_M_up + LDOS_M_down) 
plt.plot(omega*27211.6, LDOS_M_up, label = 'up')  
plt.plot(omega*27211.6, LDOS_M_down, label = 'down')
plt.title('Site %i' %M)     
plt.legend()    
    
plt.figure(3)
plt.title('Component 1')
plt.plot(u_1**2, '.', label = 'u')
plt.plot(v_1**2, '.', label = 'v')
plt.legend()

plt.figure(4)
plt.title('Component %i' %M)
plt.plot(range(len(u_M)), u_M**2, label = 'u')
plt.plot(range(len(v_M)), v_M**2, label = 'v')
plt.legend()


t2 = time.time()

print('Program finished after', (t2 - t1)/60, 'mins')



