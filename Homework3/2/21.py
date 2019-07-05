'''
Ex2.1 Homework3
function definition file
Author = Alejandro Rios
AeroStructural Optimization
ITA
'''
################################################################################
# IMPORT MODULES
import numpy as np
from llt_module import llt_module as llt
from llt_module_d import llt_module_d as llt_d
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import root
################################################################################
# EXECUTION

# Define number of vortex
n_vort = 4

# Define kynks position
x = np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 2.0, 4.0, 6.0, 8.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')
# Define gama, alpha0, chords, cl0 local and cla local
gama = np.array([1.0, 2.0, 2.0, 1.0])
alpha0 = np.array([0.0, 0.0, 0.0, 0.0])
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.array([0.0, 0.0, 0.0, 0.0])
cla = np.array([6.283, 6.283, 6.283, 6.283])

# Define ângle of attack
alpha = 5.498*np.pi/180
# Define air desnity
rho = 1
# Define velocity component vector
vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
################################################################################
# Define function that computes residuals for the given inputs
res, Sref, CL, CD, L, D = llt.tapenade_main(x, gama, alpha0, chords, cl0, cla, vinf, rho)
print('Residuals: ',res)
print('Sref: ',Sref)
print('CL: ',CL)
print('CD: ',CD)
print('L: ',L)
print('D: ',D)
