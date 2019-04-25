#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 10:53:57 2019

@author: cristina
"""
import numpy as np
import scipy.spatial.distance

N=3
M=3
T=2
a_interatomic = 1

#######calculo de distancias
i = np.arange(N+T+M)
j = np.arange(1)
I, J= np.meshgrid(j, i)
ii = np.reshape(I, ((N+T+M), ))
jj = np.reshape(J, ((N+T+M), ))
ij = zip(jj,ii)
ij = list(ij)
IJ = np.array(ij, dtype = 'double')
rr = scipy.spatial.distance.cdist(IJ, IJ, metric='euclidean')*a_interatomic#distance between sites
rr[np.where(rr == 0)] = 100   # avoid 1 / 0 errors !!

g = np.zeros([(N + M + T),(N + M + T), 4, 4], dtype=complex)

factor_diag = 7
for i in range(N):
    g_i = i
    for j in range(N):
        g_j = j
        factor = - rr[g_i, g_j]
                    
        g[g_i, g_j, 0, 0] =  factor
        g[g_i, g_j, 1, 1] =  factor
        g[g_i, g_j, 2, 2] =  factor
        g[g_i, g_j, 3, 3] =  factor
                    
        g[g_i, g_j, 0, 3] = factor
        g[g_i, g_j, 1, 2] = factor
        g[g_i, g_j, 2, 1] = factor
        g[g_i, g_j, 3, 0] = factor
        
for g_i in range(N):
        
            
    g[g_i, g_i, 0, 0] = factor_diag
    g[g_i, g_i, 1, 1] = factor_diag
    g[g_i, g_i, 2, 2] = factor_diag
    g[g_i, g_i, 3, 3] = factor_diag
    g[g_i, g_i, 0, 3] = -factor_diag
    g[g_i, g_i, 1, 2] = -factor_diag
    g[g_i, g_i, 2, 1] = -factor_diag
    g[g_i, g_i, 3, 0] = -factor_diag



factor_diag = 7
for i in range(N+T,N+T+M):
    g_i = i
    for j in range(N+T,N+T+M):
        g_j = j
        factor = - rr[g_i, g_j]
                    
        g[g_i, g_j, 0, 0] =  factor
        g[g_i, g_j, 1, 1] =  factor
        g[g_i, g_j, 2, 2] =  factor
        g[g_i, g_j, 3, 3] =  factor
                    
        g[g_i, g_j, 0, 3] = factor
        g[g_i, g_j, 1, 2] = factor
        g[g_i, g_j, 2, 1] = factor
        g[g_i, g_j, 3, 0] = factor
        
for g_i in range(N+T,N+T+M):
        
            
    g[g_i, g_i, 0, 0] = factor_diag
    g[g_i, g_i, 1, 1] = factor_diag
    g[g_i, g_i, 2, 2] = factor_diag
    g[g_i, g_i, 3, 3] = factor_diag
    g[g_i, g_i, 0, 3] = -factor_diag
    g[g_i, g_i, 1, 2] = -factor_diag
    g[g_i, g_i, 2, 1] = -factor_diag
    g[g_i, g_i, 3, 0] = -factor_diag
    
    
    
    
    
    
    
    
    
    