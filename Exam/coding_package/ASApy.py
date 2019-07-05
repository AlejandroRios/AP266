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
from asa_module import asa_module  as asamod
from scipy.optimize import root
################################################################################
# EXECUTION

# Define number of vortex
n_panels = 4

#################################################################################
# GEOMETRIC PARAMETERS

# Define kynks position
Xfem = np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 2.0, 4.0, 6.0, 8.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')

#################################################################################
# AED PARAMETERS
# Define gama, alpha0, chords, cl0 local and cla local
gama = np.array([80.0, 80.0, 80.0, 80.0])
alpha0 = np.array([0.0, 0.0, 0.0, 0.0])
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.array([0.0, 0.0, 0.0, 0.0])
cla = np.array([6.283, 6.283, 6.283, 6.283])
# Velocity modulos
Vinfm = 60.0
# Define ângle of attack
alpha = 5*np.pi/180
# Define air desnity
rhoAir = 1.225
# Define velocity component vector
Vinf = np.array([Vinfm*np.cos(alpha), 0.0, Vinfm*np.sin(alpha)])

#################################################################################
# STUC PARAMETERS

R = np.array([0.1, 0.1, 0.1, 0.1])
t = np.array([0.005, 0.005, 0.005, 0.005])
E = 73.1e9
rhoMat = 2780.0
sigmaY = 324.0e6
pKS = 200.0
conIDs = np.array([5,6], dtype = 'int32')
d = np.array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])

#################################################################################
# ADDITIONAL PARAMETERS

CD0 = 0.0270
fixedMass = 700.0
g = 9.8
Endurance = 4.0*60.0*60.0
TSFC = 0.5/3600.0
loadFactor = 3.0*1.5

################################################################################
# Define function that computes residuals for the given inputs

resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL= asamod.asa_analysis(gama, alpha0, chords, Xfem, R, t, d, cl0, cla, Vinf, rhoAir, E, rhoMat, sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor)


print('Residuals llt: ',resllt)
print('Residuals fem: ',resfem)
print('lift Exc.: ',liftExcess)
print('margins: ',margins)
print('KSmargin: ',KSmargin)
print('FB: ', FB)
print('Weight: ',Weight)
print('Sref: ',Sref)
print('CL: ',CL)




