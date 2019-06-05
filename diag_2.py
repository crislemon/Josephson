#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:17:14 2019

@author: cristina
"""

import numpy as np
from numpy import linalg as LA
diag = LA.eigh
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
import time

pi = np.pi
#import pybinding as pb

t1 = time.time()

N = 2000 #number of sites
m = 1.0 #effective mass
delta =1.35/27211.6 #SC gap
mu = 1.0/27211.6 #chemical potential
mu = 0.0
a = 4.98/0.529 ##lattice constant
#e0 = 2.0/27211.6
#t = 0.0

H = np.zeros([N * 2,N  * 2], dtype=complex)
h = np.zeros([N, N, 2, 2], dtype=complex)

factor = 1/(m*a**2) - mu
factor_2 = -1/(2*m*a**2)

for i in range(N):
    g_i = i
    
    #diagonal on atom index                
    h[g_i, g_i, 0, 1] = delta
    h[g_i, g_i, 1, 0] = delta    
                    
    h[g_i, g_i, 0, 0] = factor
    h[g_i, g_i, 1, 1] = - factor
   
for i in range(N-1):
    g_i = i
    g_j = i + 1    
    
    #off-diagonal terms    
    h[g_i, g_j, 0, 0] = factor_2
    h[g_i, g_j, 1, 1] = - factor_2   
    
    h[g_j, g_i, 0, 0] = factor_2
    h[g_j, g_i, 1, 1] = - factor_2
    

for i in range(N):
        for j in range(N):
            for t_i in range(2):
                for t_j in range(2):
                    H[(i) * 2 + t_i, (j) * 2 + t_j] = h[i, j, t_i, t_j]
                    
H = np.matrix(H)  
T = np.allclose(H, H.getH())###check if Hermitian
print(T)

(E, psi) = diag(H)####diagonalize H
u_1 = np.zeros(len(E))
v_1 = np.zeros(len(E))

u_M = np.zeros(len(E))
v_M = np.zeros(len(E))


M = int(N/2.0)

for i in range(len(E)):
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
    

###calculate LDOS
omega = np.linspace(-4*delta, 4*delta, 2000)#omega vector     

LDOS_1_up = np.zeros(len(omega))
LDOS_1_down = np.zeros(len(omega))

LDOS_M_up = np.zeros(len(omega))
LDOS_M_down = np.zeros(len(omega))

D = 0.08/27211.6
for i in range(len(omega)):
    LDOS_1_up[i] = LDOS_up(omega[i], E, u_1, D)    
    LDOS_1_down[i] = LDOS_down(omega[i], E, v_1, D)

    LDOS_M_up[i] = LDOS_up(omega[i], E, u_M, D)    
    LDOS_M_down[i] = LDOS_down(omega[i], E, v_M, D)
    
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
plt.plot(u_M**2, '.', label = 'u')
plt.plot(v_M**2, '.', label = 'v')
plt.legend()


t2 = time.time()

print('Program finished after', t2 - t1)



