#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:08:06 2019

@author: cristina
"""

import numpy as np
from numpy import linalg as LA
diag = LA.eigh
import matplotlib.pyplot as plt
import time
plt.rcParams.update({'font.size': 13})

t1 = time.time()

pi = np.pi
N = 3000 #number of sites
m = 1.0 #effective mass
delta =1.35/27211.6 #SC gap
mu = 1.0/27211.6 #chemical potential
mu = 0.0
a = 4.98/0.529 ##lattice constant

factor = 1/(m*a**2) - mu
factor_2 = -1/(2*m*a**2)


H = np.zeros([N * 4,N  * 4], dtype=complex)
h = np.zeros([N, N, 4, 4], dtype=complex)

for i in range(N):
    g_i = i
    g_j = i + 1
    
    
    #diagonal on atom index                
    h[g_i, g_i, 0, 3] = -delta
    h[g_i, g_i, 1, 2] = delta
    h[g_i, g_i, 2, 1] = delta 
    h[g_i, g_i, 3, 0] = -delta
                    
    h[g_i, g_i, 0, 0] = factor
    h[g_i, g_i, 1, 1] = factor
    h[g_i, g_i, 2, 2] = - factor
    h[g_i, g_i, 3, 3] = - factor
    
for i in range(N-1):
    g_i = i
    g_j = i + 1

    #off-diagonal terms
    
    h[g_i, g_j, 0, 0] = factor_2
    h[g_i, g_j, 1, 1] = factor_2
    h[g_i, g_j, 2, 2] = - factor_2
    h[g_i, g_j, 3, 3] = - factor_2
    
    h[g_j, g_i, 0, 0] = factor_2
    h[g_j, g_i, 1, 1] = factor_2
    h[g_j, g_i, 2, 2] = - factor_2
    h[g_j, g_i, 3, 3] = - factor_2
    

for i in range(N):
        for j in range(N):
            for t_i in range(4):
                for t_j in range(4):

                    H[(i) * 4 + t_i, (j) * 4 + t_j] = h[i, j, t_i, t_j]
H = np.matrix(H)  
T = np.allclose(H, H.getH())###check if Hermitian
print('Is H Hermitian matrix?', T)
(E, psi) = diag(H)####diagonalize H

u_1up = np.zeros(len(E))
u_1down = np.zeros(len(E))
v_1up = np.zeros(len(E))
v_1down = np.zeros(len(E))

u_Mup = np.zeros(len(E))
u_Mdown = np.zeros(len(E))
v_Mup = np.zeros(len(E))
v_Mdown = np.zeros(len(E))


M = int(N/2.0)

for i in range(len(E)):
    u_1up[i] = psi[0,i]
    u_1down[i] = psi[1,i]
    v_1up[i] = psi[2,i]
    v_1down[i] = psi[3,i]
    
    u_Mup[i] = psi[4*(M-1),i]
    u_Mdown[i] = psi[4*M - 3,i]
    v_Mup[i] = psi[4*M - 2,i]
    v_Mdown[i] = psi[4*M - 1,i]

####calculate DOS
def LDOS_up(omega, E, u, Damping):
    t = sum ( u**2 / (omega - E + 1j*Damping) )
    tt = -1/pi*np.imag(t)
    return(tt)
    
def LDOS_down(omega, E, v, Damping):
    t = sum ( v**2 / (omega + E + 1j*Damping) )
    tt = -1/pi*np.imag(t)
    return(tt)

omega = np.linspace(-4*delta, 4*delta, 2000)     
    
LDOS_1_u_up = np.zeros(len(omega))
LDOS_1_u_down = np.zeros(len(omega))
LDOS_1_v_up = np.zeros(len(omega))
LDOS_1_v_down = np.zeros(len(omega))

LDOS_M_u_up = np.zeros(len(omega))
LDOS_M_u_down = np.zeros(len(omega))
LDOS_M_v_up = np.zeros(len(omega))
LDOS_M_v_down = np.zeros(len(omega))

D = 0.08/27211.6
for i in range(len(omega)):
    LDOS_1_u_up[i] = LDOS_up(omega[i], E, u_1up, D)    
    LDOS_1_u_down[i] = LDOS_up(omega[i], E, u_1down, D)
    LDOS_1_v_up[i] = LDOS_up(omega[i], E, v_1up, D)    
    LDOS_1_v_down[i] = LDOS_up(omega[i], E, v_1down, D)

    LDOS_M_u_up[i] = LDOS_up(omega[i], E, u_Mup, D)    
    LDOS_M_u_down[i] = LDOS_up(omega[i], E, u_Mdown, D)
    LDOS_M_v_up[i] = LDOS_up(omega[i], E, v_Mup, D)    
    LDOS_M_v_down[i] = LDOS_up(omega[i], E, v_Mdown, D)
    
plt.figure(1)
plt.plot(omega*27211.6, LDOS_1_u_up + LDOS_1_v_down, label = 'u up + v down') 
plt.plot(omega*27211.6, LDOS_1_u_up, label = 'u up')
plt.plot(omega*27211.6, LDOS_1_v_down, label = 'v down')
plt.title('Site 1')     
plt.legend()

plt.figure(2)
plt.plot(omega*27211.6, LDOS_1_u_down + LDOS_1_v_up, label = 'u down + v up')   
plt.plot(omega*27211.6, LDOS_1_u_down, 'r', label = 'u down')
plt.plot(omega*27211.6, LDOS_1_v_up, 'g', label = 'v up')  

plt.title('Site 1')     
plt.legend()        
    
plt.figure(3)
plt.plot(omega*27211.6, LDOS_M_u_up + LDOS_M_v_down, label = 'u up + v down') 
plt.plot(omega*27211.6, LDOS_M_u_up, label = 'u up')    
plt.plot(omega*27211.6, LDOS_M_v_down, label = 'v down')
plt.title('Site %i' %M)     
plt.legend()

plt.figure(4)
plt.plot(omega*27211.6, LDOS_M_u_down + LDOS_M_v_up, label = 'u down + v up') 
plt.plot(omega*27211.6, LDOS_M_u_down, 'r', label = 'u down')
plt.plot(omega*27211.6, LDOS_M_v_up, 'g', label = 'v up')  
plt.title('Site %i' %M)     
plt.legend()        
    
plt.figure(5)
plt.title('Component 1')
plt.plot(u_1up**2, '.', label = 'u up')
plt.plot(v_1down**2, '.', label = 'v down')
plt.legend()

plt.figure(6)
plt.title('Component 1')
plt.plot(u_1down**2, 'r.', label = 'u down')
plt.plot(v_1up**2, 'g.', label = 'v up')
plt.legend()

plt.figure(7)
plt.title('Component %i' %M)
plt.plot(u_Mup**2, '.', label = 'u up')
plt.plot(v_Mdown**2, '.', label = 'v down')
plt.legend()

plt.figure(8)
plt.title('Component %i' %M)
plt.plot(u_Mdown**2, 'r.', label = 'u down')
plt.plot(v_Mup**2, 'g.', label = 'v up')
plt.legend()


t2 = time.time()

print('Program finished after', (t2 - t1)/60, 'min')












