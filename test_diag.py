#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:37:05 2019

@author: cristina
"""

import numpy as np
from itertools import chain

N = 3
M = 2
h = np.zeros([2*N + M, 2*N + M])
range1 = range(N)
range2 = range(N+M, 2*N+M)
range_total = chain(range1, range2)

factor = 3.0
factor_2 = 2.0
factor_3 = 5.0

#diagonal terms
for i in range_total:
    g_i = i
    
    h[g_i, g_i] = factor
    h[g_i, g_i] = factor
   

range1_offdiagonal = range(N - 1)
range2_offdiagonal = range(N+M, 2*N+M - 1)
range_offdiagonal = chain(range1_offdiagonal, range2_offdiagonal)

#offdiagonal terms
for i in range_offdiagonal:
    g_i = i
    g_j = i + 1
    
    h[g_i, g_j] = factor_2
    h[g_j, g_i] = - factor_2
    
    
#hopping between the 2 Chains
h[N - 1, N + M] = factor_3
h[N + M, N - 1] = - factor_3