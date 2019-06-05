#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 16:03:44 2018

@author: cristina
"""

import numpy as np
import scipy.spatial.distance

sin = np.sin
cos = np.cos
exp =np.exp
ui = complex(0.0, 1.0)

#hoppig between electrodes


def Self_Energy(N, M, T, a_interatomic, t, phi1, phi2):
    
     Self = np.zeros([N + T + M, N + T + M, 4, 4], dtype=complex)
     Self2 = np.zeros([(N + T + M) * 4, (N + T + M) * 4], dtype=complex)
     
     t = t*exp(1j*(phi1-phi2))
     
     i = np.arange(N+T+M)
     j = np.arange(1)
     I, J= np.meshgrid(j, i)
     ii = np.reshape(I, ((N+T+M), ))
     jj = np.reshape(J, ((N+T+M), ))
     ij = zip(jj,ii)
     ij = list(ij)
     IJ = np.array(ij, dtype = 'double')
     rr = scipy.spatial.distance.cdist(IJ, IJ, metric='euclidean')*a_interatomic#distance between sites

     
     '''hopping between the two electrodes'''
     g_i = N - 1
     g_j = N + T
     
#     Self [g_i, g_j, 0, 0]= t/rr[g_i,g_j]
#     Self [g_i, g_j, 1, 1]= t/rr[g_i,g_j]
#     Self [g_i, g_j, 2, 2]= t/rr[g_i,g_j]
#     Self [g_i, g_j, 3, 3]= t/rr[g_i,g_j]
#     
#     Self [g_j, g_i, 0, 0]= t/rr[g_j, g_i]
#     Self [g_j, g_i, 1, 1]= t/rr[g_j, g_i]
#     Self [g_j, g_i, 2, 2]= t/rr[g_j, g_i]
#     Self [g_j, g_i, 3, 3]= t/rr[g_j, g_i]
     
     Self [g_i, g_j, 0, 2]= t/rr[g_i,g_j]
     Self [g_i, g_j, 1, 3]= t/rr[g_i,g_j]
     Self [g_i, g_j, 2, 0]= t/rr[g_i,g_j]
     Self [g_i, g_j, 3, 1]= t/rr[g_i,g_j]
     
     Self [g_j, g_i, 0, 2]= t/rr[g_j, g_i]
     Self [g_j, g_i, 1, 3]= t/rr[g_j, g_i]
     Self [g_j, g_i, 2, 0]= t/rr[g_j, g_i]
     Self [g_j, g_i, 3, 1]= t/rr[g_j, g_i]
     
            
        
     
            
     for i in range((N+T+M)):
        for j in range((N+T+M)):
            for t_i in range(4):
                for t_j in range(4):
                    Self2[(i) * 4 + t_i, (j) * 4 + t_j] = Self[i, j, t_i, t_j]
                    
                    
     return(Self2)
            