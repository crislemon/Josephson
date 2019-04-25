#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 17:53:46 2018

@author: cristina
"""

"this program makes all the plots"


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import time
import matplotlib.pyplot as plt

#parameters
pi=np.pi
d = 1.0 #distance between sites
k_F = 0.183
DOS = 1.0

delta = 0.75/27211.6 #SC gap
N_omega = 2003
range_omega = 4


N = 7
M = 7
T = 3

t= 1.5

################################################# We solve Dyson's equation

import Green_Chain as GC
t1=time.time()
(gg , N, M, T, N_omega , vv, Go, Self2) = GC.Josephson(d, k_F, DOS, delta, N_omega, range_omega, t, N, M, T)
t2 = time.time()
 
print('The program is finished after', t2 - t1)


##################################################


"The spectrum is obtained all Nambu components"
spectro = np.zeros([(N + T + M), N_omega], dtype= 'float')


for i in range((N + T + M)):
         #I = i_atom*N_x + j_atom
         for i_omega in range(N_omega):
             
             
             tr3 = gg[i*4 + 0, i*4 + 0, i_omega] + gg[i*4 + 1, i*4 + 1, i_omega]
             spectro[i, i_omega] = - (tr3.imag)/(pi)
             
             
             
#####
"Plot the spectrum in the border of both SC"

spectro_1 = spectro[N-1, :]#spectrum in the border
spectro_2 = spectro[N+T, :]#spectrum in the border
spectro_3 = spectro[N-2, :]#spectrum in "bulk"
plt.figure(1)
plt.plot(vv, spectro_1, linewidth=1.0, label = 'total')
plt.xlabel('meV')
plt.ylabel('PDOS')


plt.figure(2)
plt.plot(vv, spectro_2, linewidth=1.0, label = 'total')
plt.xlabel('meV')
plt.ylabel('PDOS')

plt.figure(3)
plt.plot(vv, spectro_3, linewidth=1.0, label = 'total')
plt.xlabel('meV')
plt.ylabel('PDOS')







