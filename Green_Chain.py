#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 12:37:27 2018

@author: cristina
"""
#everything in atomic units
import numpy as np

def Josephson(nstep, k_f, DOS, delta, N_omega, range_omega, t, N, M, T, phi_1, phi_2):
    
    # d = nstep*a distance between sites
    
    "SC electrodes" 
    #N = number of sites in first electrode
    #M = number of sites in second electrode
    #T = number of sites in vacuum
    
    "Material data Bi2Pd"
    Damping = 0.02/27211.6 #Dynes damping
    Delta=0.75/27211.6 #SC gap
    Delta = delta
    DOS_o = DOS #Normal phase DOS
    Fermi_k = k_f
    mass_eff=1 #SC Band effective mass
    a_interatomic=nstep*3.36/0.529

    

    "we define the omega vector"
    #from -N_delta/2 to N_delta/2
    N_delta = range_omega
    
    Romega = np.zeros([N_omega])
    Romega=np.array(Romega, np.longdouble)
    step_omega=N_delta*Delta/(N_omega-1)

    for i_omega in range(N_omega):
        Romega[i_omega] = (-N_delta/2.*Delta+(i_omega)*step_omega)
     
        
    Romega = np.array(Romega)
    vv=Romega*27211.6

    
    "We calculate the Green's functions and solve Dyson eq"
    
    
    import Self_Energy_loop as SL
    Self2 = SL.Self_Energy(N, M, T, a_interatomic, t)
    
    GG = np.zeros([4 * (N + T + M) , 4 * (N + T + M), N_omega], dtype=complex)
    
    for i_omega in range(N_omega):
        
        omega = Romega[i_omega]
    
         #BCS Green's function
        import Free_Gree_loop as FG
        Go = FG.Free_Green(N, M, T, omega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic, phi_1, phi_2)
        
        #Solve Dyson's equation
        import Dyson as Dy
        gg = Dy.Dyson_eq(Go, Self2, N, M, T)
        
        GG[:,:, i_omega] = gg
        
        
    return(GG, N, M, T, N_omega , vv, Go, Self2)    
        
        
        
    

