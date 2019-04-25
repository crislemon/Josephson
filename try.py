#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:50:30 2019

@author: cristina
"""


d =1.0
k_F = 0.183
U=0
j=1.0
DOS=1.0
s=3
delta = 0.75/27211.6 #SC gap
N_omega = 2003
range_omega = 4
import Josephson as JS
(N, M, T, N_omega , vv, Go, Self2) = JS.Josephson(d, k_F, U, j, DOS, s, delta, N_omega, range_omega)

