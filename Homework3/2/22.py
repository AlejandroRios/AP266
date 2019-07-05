'''
Ex2.2 Homework3
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
# Define input seed for reversal AD
resb = np.array([5.0, -0.3, 0.0, 1.0])
srefb = -0.3
clb = 1.3
cdb = -0.5
lb = 3.0
db = 0.02

# Initializing input seedses)
xb =np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')

gamab = np.array([0.0, 0.0, 0.0, 0.0])
alpha0b = np.array([0.0, 0.0, 0.0, 0.0])
chordsb = np.array([0.0, 0.0, 0.0, 0.0])
gama0 = np.array([0.0, 0.0, 0.0, 0.0])

# Define function that computes residuals for given gama
def resfunc(gama) :
 res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
 return res

# Solve the AED pŕoblem
sol = root(resfunc,gama0)

# Get the gama values that solve LLT problem
gama = sol.x
print('Gama:',gama)
# Get residuals for solved gama
res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)

# Compute functions for the gama that solves LLT problem
Sref, CL, CD, L, D = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)

print('Sref ',Sref)
print('CL: ',CL)
print('CD: ',CD)
print('L: ',L)
print('D: ',D)
