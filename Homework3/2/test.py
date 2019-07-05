'''
03_solve_adjoint.py
This script solves the LLT problem and solves the adjoint method
'''

# IMPORTS
import numpy as np
from llt_module import llt_module as llt
from llt_module_d import llt_module_d as llt_d
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import minimize
from scipy.optimize import root
import matplotlib.pyplot as plt

# Define number of bars
n_vort = 10
b = 8
chord = 1
cla = 6.283

x = np.zeros([3,n_vort])
x[0] = x[0,:]
x[1] =  np.linspace(0, b, n_vort)
x[2] = x[2,:]
x = np.array(x,order='F')

gama = np.ones([1,n_vort])*2
gama[0,0] = 1
gama[-1,-1] = 1


alpha0 = np.zeros([1,n_vort])

chords = np.ones([1,n_vort])*chord

cl0 = np.zeros([1,n_vort])

cla = np.ones([1,n_vort])*cla
