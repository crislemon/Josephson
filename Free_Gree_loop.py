#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:32:08 2018

@author: cristina
"""

import numpy as np
import scipy.spatial.distance

# Functions
pi = np.pi
sin = np.sin
cos = np.cos
sqrt = np.sqrt
exp = np.exp


def Free_Green(N, M, T, lomega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic):

    G = np.zeros([(N + M + T) * 4,(N + M + T) * 4], dtype=complex)
    g = np.zeros([(N + M + T),(N + M + T), 4, 4], dtype=complex)
    
    omega = lomega + 1j * Damping

    
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
    
    
    '''First SC electrode'''
    # Non diagonal in atom
    SS = sqrt(Delta**2 - omega**2)
    xi = Fermi_k / (mass_eff * SS)
    
    
    for i in range(N):
            g_i = i
            for j in range(N):
                g_j = j
                factor = - pi * DOS_o * exp(-rr[g_i, g_j]/ xi) / (SS * Fermi_k * rr[g_i, g_j])
                    
                g[g_i, g_j, 0, 0] = ( omega * sin(Fermi_k * rr[g_i, g_j]) + SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                g[g_i, g_j, 1, 1] = ( omega * sin(Fermi_k * rr[g_i, g_j]) - SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                g[g_i, g_j, 2, 2] = ( omega * sin(Fermi_k * rr[g_i, g_j]) + SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                g[g_i, g_j, 3, 3] = ( omega * sin(Fermi_k * rr[g_i, g_j]) - SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                    
                g[g_i, g_j, 0, 3] = - Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
                g[g_i, g_j, 1, 2] = Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
                g[g_i, g_j, 2, 1] = Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
                g[g_i, g_j, 3, 0] = - Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
    

    # Diagonal in atom
    omega = lomega + 1j * Damping
    SS = sqrt(Delta**2 - omega**2)
    factor_diag = - pi * DOS_o / SS

    for g_i in range(N):
        
            
            g[g_i, g_i, 0, 0] = omega * factor_diag
            g[g_i, g_i, 1, 1] = omega * factor_diag
            g[g_i, g_i, 2, 2] = omega * factor_diag
            g[g_i, g_i, 3, 3] = omega * factor_diag
            g[g_i, g_i, 0, 3] = -Delta * factor_diag
            g[g_i, g_i, 1, 2] = Delta * factor_diag
            g[g_i, g_i, 2, 1] = Delta * factor_diag
            g[g_i, g_i, 3, 0] = -Delta * factor_diag



    '''Second SC electrode'''
    # Non diagonal in atom
    SS = sqrt(Delta**2 - omega**2)
    xi = Fermi_k / (mass_eff * SS)
    
    
    for i in range(N+T,N+T+M):
        g_i = i
        for j in range(N+T,N+T+M):
            g_j = j
            factor = - pi * DOS_o * exp(-rr[g_i, g_j]/ xi) / (SS * Fermi_k * rr[g_i, g_j])
                    
            g[g_i, g_j, 0, 0] = ( omega * sin(Fermi_k * rr[g_i, g_j]) + SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
            g[g_i, g_j, 1, 1] = ( omega * sin(Fermi_k * rr[g_i, g_j]) - SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
            g[g_i, g_j, 2, 2] = ( omega * sin(Fermi_k * rr[g_i, g_j]) + SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
            g[g_i, g_j, 3, 3] = ( omega * sin(Fermi_k * rr[g_i, g_j]) - SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                    
            g[g_i, g_j, 0, 3] = - Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
            g[g_i, g_j, 1, 2] = Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
            g[g_i, g_j, 2, 1] = Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
            g[g_i, g_j, 3, 0] = - Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
    

    # Diagonal in atom
    omega = lomega + 1j * Damping
    SS = sqrt(Delta**2 - omega**2)
    factor_diag = - pi * DOS_o / SS

    for g_i in range(N+T,N+T+M):
        
            
            g[g_i, g_i, 0, 0] = omega * factor_diag
            g[g_i, g_i, 1, 1] = omega * factor_diag
            g[g_i, g_i, 2, 2] = omega * factor_diag
            g[g_i, g_i, 3, 3] = omega * factor_diag
            g[g_i, g_i, 0, 3] = -Delta * factor_diag
            g[g_i, g_i, 1, 2] = Delta * factor_diag
            g[g_i, g_i, 2, 1] = Delta * factor_diag
            g[g_i, g_i, 3, 0] = -Delta * factor_diag
            

    for i in range(N+M+T):
        for j in range(N+M+T):
            for t_i in range(4):
                for t_j in range(4):
                    G[(i) * 4 + t_i, (j) * 4 + t_j] = g[i, j, t_i, t_j]
    
    
    
    return (G)
