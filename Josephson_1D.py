#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:04:05 2019

@author: cristina
"""

import numpy as np
from itertools import chain
from numpy import linalg as LA
diag = LA.eigh
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})
import time

pi = np.pi
exp = np.exp
t1 = time.time()

N = 2000 #number of sites
M = 200 #number of empty sites
m = 1.0 #effective mass
delta =1.35/27211.6 #SC gap
mu = 1.0/27211.6 #chemical potential
mu = 0.0
a = 4.98/0.529 ##lattice constant
phi = pi/2.0#phase of second SC
phi = 0.0

H = np.zeros([2*(2*N + M), 2*(2*N + M)], dtype=complex)
h = np.zeros([2*N + M, 2*N + M, 2, 2], dtype=complex)

factor = 1/(m*a**2) - mu
factor_2 = -1/(2*m*a**2)
hopping = factor_2*10
hopping = 0.0

#diagonal terms
range1_diagonal = range(N)
range2_diagonal = range(N+M, 2*N+M - 1)

for i in range1_diagonal:
    g_i = i
    
    h[g_i, g_i, 0, 1] = delta
    h[g_i, g_i, 1, 0] = delta    
                    
    h[g_i, g_i, 0, 0] = factor
    h[g_i, g_i, 1, 1] = - factor
    
    
for i in range2_diagonal:
    g_i = i
    
    h[g_i, g_i, 0, 1] = delta*exp(1j*phi)
    h[g_i, g_i, 1, 0] = delta*exp(-1j*phi)    
                    
    h[g_i, g_i, 0, 0] = factor
    h[g_i, g_i, 1, 1] = - factor

#off - diagonal terms
range1_offdiagonal = range(N - 1)
range2_offdiagonal = range(N+M, 2*N+M - 1)
range_offdiagonal = chain(range1_offdiagonal, range2_offdiagonal)

for i in range_offdiagonal:
    g_i = i
    g_j = i + 1
    
    h[g_i, g_j, 0, 0] = factor_2
    h[g_i, g_j, 1, 1] = - factor_2   
    
    h[g_j, g_i, 0, 0] = factor_2
    h[g_j, g_i, 1, 1] = - factor_2
    

#hopping between the 2 Chains
h[N - 1, N + M, 0, 0] = hopping
h[N - 1, N + M, 1, 1] = - hopping

h[N + M, N - 1, 0, 0] = hopping
h[N + M, N - 1, 1, 1] = - hopping

for i in range(2*N + M):
        for j in range(2*N + M):
            for t_i in range(2):
                for t_j in range(2):
                    H[(i) * 2 + t_i, (j) * 2 + t_j] = h[i, j, t_i, t_j]
                    
H = np.matrix(H)  
T = np.allclose(H, H.getH())###check if Hermitian
print('Is H an Hermitian matrix?', T)

(E, psi) = diag(H)####diagonalize H


####LDOS functions
def LDOS_up(omega, E, u, Damping):
    t = sum ( u**2 / (omega - E + 1j*Damping) )
    tt = -1/pi*np.imag(t)
    return(tt)
    
def LDOS_down(omega, E, v, Damping):
    t = sum ( v**2 / (omega + E + 1j*Damping) )
    tt = -1/pi*np.imag(t)
    return(tt)



#### u and v components in the Nth atom
u_borde1 = np.zeros(len(E))
v_borde1 = np.zeros(len(E))
I = N - 1

u_borde2 = np.zeros(len(E))
v_borde2 = np.zeros(len(E))
I2 = N + M - 1

u_bulk1 = np.zeros(len(E))
v_bulk1 = np.zeros(len(E))
I3 = int(N/2) - 1

u_bulk2 = np.zeros(len(E))
v_bulk2 = np.zeros(len(E))
I4 = N + M + int(N/2.0) - 1

I = N 
for i in range(len(E)):
    u_borde1[i] = psi[2*I-2,i]
    v_borde1[i] = psi[2*I-1,i]
    
    u_borde2[i] = psi[2*I2-2,i]
    v_borde2[i] = psi[2*I2-1,i]
    
    u_bulk1[i] = psi[2*I3-2,i]
    v_bulk1[i] = psi[2*I3-1,i]
    
    u_bulk2[i] = psi[2*I4-2,i]
    v_bulk2[i] = psi[2*I4-1,i]

###calculate LDOS
omega = np.linspace(-4*delta, 4*delta, 2000)#omega vector     

LDOS_borde1_up = np.zeros(len(omega))
LDOS_borde1_down = np.zeros(len(omega))

LDOS_borde2_up = np.zeros(len(omega))
LDOS_borde2_down = np.zeros(len(omega))

LDOS_bulk1_up = np.zeros(len(omega))
LDOS_bulk1_down = np.zeros(len(omega))

LDOS_bulk2_up = np.zeros(len(omega))
LDOS_bulk2_down = np.zeros(len(omega))

D = 0.02/27211.6
for i in range(len(omega)):

    LDOS_borde1_up[i] = LDOS_up(omega[i], E, u_borde1, D)    
    LDOS_borde1_down[i] = LDOS_up(omega[i], E, v_borde1, D)
    
    LDOS_borde2_up[i] = LDOS_up(omega[i], E, u_borde2, D)    
    LDOS_borde2_down[i] = LDOS_up(omega[i], E, v_borde2, D)
    
    LDOS_bulk1_up[i] = LDOS_up(omega[i], E, u_bulk1, D)    
    LDOS_bulk1_down[i] = LDOS_up(omega[i], E, v_bulk1, D)
    
    LDOS_bulk2_up[i] = LDOS_up(omega[i], E, u_bulk2, D)    
    LDOS_bulk2_down[i] = LDOS_up(omega[i], E, v_bulk2, D)


###plot LDOS    
plt.figure(1)
plt.plot(omega*27211.6, LDOS_borde1_up + LDOS_borde1_down) 
plt.plot(omega*27211.6, LDOS_borde1_up, label = 'up')  
plt.plot(omega*27211.6, LDOS_borde1_down, label = 'down')
plt.title('Borde SC 1')
#plt.title('Site %i' %I)     
plt.legend()   

plt.figure(2)
plt.plot(omega*27211.6, LDOS_borde2_up + LDOS_borde2_down) 
plt.plot(omega*27211.6, LDOS_borde2_up, label = 'up')  
plt.plot(omega*27211.6, LDOS_borde2_down, label = 'down')
plt.title('Borde SC 2')
#plt.title('Site %i' %I)     
plt.legend()   

plt.figure(3)
plt.plot(omega*27211.6, LDOS_bulk1_up + LDOS_bulk1_down) 
plt.plot(omega*27211.6, LDOS_bulk1_up, label = 'up')  
plt.plot(omega*27211.6, LDOS_bulk1_down, label = 'down')
plt.title('Bulk SC 1')
#plt.title('Site %i' %I)     
plt.legend()

plt.figure(4)
plt.plot(omega*27211.6, LDOS_bulk2_up + LDOS_bulk2_down) 
plt.plot(omega*27211.6, LDOS_bulk2_up, label = 'up')  
plt.plot(omega*27211.6, LDOS_bulk2_down, label = 'down')
plt.title('Bulk SC 2')
#plt.title('Site %i' %I)     
plt.legend()    



t2 = time.time()
print('Program finished after', (t2 - t1)/60.0, 'mins')








