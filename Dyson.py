# -*- coding: utf-8 -*-

import numpy as np
from numpy.linalg import inv

#we solve Dyson's equation
def Dyson_eq(Go, Self2, N, M, T):
    
    
    Id = np.identity(4 * (N + M + T))
   
    matrx_inv = inv(Id - np.dot(Go, Self2))##### + o -?????
    gg = np.dot(matrx_inv , Go)
    
    
    
    
    
    
    return(gg)
