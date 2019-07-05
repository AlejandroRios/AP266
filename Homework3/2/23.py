'''
Ex2.3 Homework3
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
# Define input seed for reversal AD
resb = np.array([5.0, -0.3, 0.0, 1.0])
srefb = -0.3
clb = 1.3
cdb = -0.5
lb = 3.0
db = 0.02


res, sref, cl, cd, l, d = llt.tapenade_main(x, gama, alpha0, chords, cl0, cla, vinf, rho)


###############################################################################
# Define function that receives resb and computes the adjoint
# function for derivatives of CL
def adjfunbc(resb):
    # Define active input seed for reversal AD
    srefb = 0.0
    clb = 1.0
    cdb = 0.0
    lb = 0.0
    db = 0.0

    # Initializing input seedses)
    xb =np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')
    gamab = np.array([0.0, 0.0, 0.0, 0.0])
    alpha0b = np.array([0.0, 0.0, 0.0, 0.0])
    chordsb = np.array([0.0, 0.0, 0.0, 0.0])

    # Run reverse mode
    llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
    chords, chordsb, cl0, cla, vinf, rho, res, resb.copy(), sref, srefb, cl,
    clb, cd, cdb, l, lb, d, db)

    # Compute the residuals of the adjoint equation
    adj_res = gamab
    return adj_res
###############################################################################
# Solve the adjoint problem
psi0 = np.ones(n_vort)*10**6
sol = root(adjfunbc, psi0)

# Get adjoint variables (psi)
psi = sol.x

# Define active input seed for reversal AD
srefb = 0
clb = 1
cdb = 0
lb = 0
db = 0

# Initializing input seedses)
xb =np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')
gamab = np.array([0.0, 0.0, 0.0, 0.0])
alpha0b = np.array([0.0, 0.0, 0.0, 0.0])
chordsb = np.array([0.0, 0.0, 0.0, 0.0])
resb = psi.copy()

# Run reverse mode of tapenade_main
llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
clb, cd, cdb, l, lb, d, db)

# Get total derivatives for CD
dCLdX = xb
dCLdalpha0 = alpha0b
dCLdc = chordsb
rCL = gamab
print('psiCL:',psi)
print('dCLdX', dCLdX)
print('dCLdalpha0', dCLdalpha0)
print('dCLdc', dCLdc)
print('rCL', rCL)
###############################################################################
# Define function that receives resb and computes the adjoint
# function for derivatives of CD
def adjfunbc(resb):
    # Define active input seed for reversal AD
    srefb = 0.0
    clb = 0.0
    cdb = 1.0
    lb = 0.0
    db = 0.0

    # Initializing input seedses)
    xb =np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')
    gamab = np.array([0.0, 0.0, 0.0, 0.0])
    alpha0b = np.array([0.0, 0.0, 0.0, 0.0])
    chordsb = np.array([0.0, 0.0, 0.0, 0.0])

    # Run reverse mode
    llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
    chords, chordsb, cl0, cla, vinf, rho, res, resb.copy(), sref, srefb, cl,
    clb, cd, cdb, l, lb, d, db)

    # Compute the residuals of the adjoint equation
    adj_res = gamab
    return adj_res
###############################################################################
# Solve the adjoint problem
psi0 = np.ones(n_vort)*10**6
sol = root(adjfunbc, psi0)

# Get adjoint variables (psi)
psi = sol.x

# Define active input seed for reversal AD
srefb = 0
clb = 0
cdb = 1
lb = 0
db = 0

# Initializing input seeds
xb =np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')
gamab = np.array([0.0, 0.0, 0.0, 0.0])
alpha0b = np.array([0.0, 0.0, 0.0, 0.0])
chordsb = np.array([0.0, 0.0, 0.0, 0.0])
resb = psi.copy()

# Run reverse mode of tapenade_main
llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
clb, cd, cdb, l, lb, d, db)

# Get total derivatives for CD
dCDdX = xb
dCD0dalpha = alpha0b
dCDdc = chordsb
rCD = gamab
print('psiCD:',psi)
print('dCDdX', dCDdX)
print('dCDdalpha0', dCD0dalpha)
print('dCDdc', dCDdc)
print('rCD', rCD)

###############################################################################
# Define function that receives resb and computes the adjoint
# function for derivatives of Sref
# Define active input
srefb = 1.0
clb = 0.0
cdb = 0.0
lb = 0.0
db = 0.0

# Initializing input seeds
xb =np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]],order='F')

gamab = np.array([0.0, 0.0, 0.0, 0.0])
alpha0b = np.array([0.0, 0.0, 0.0, 0.0])
chordsb = np.array([0.0, 0.0, 0.0, 0.0])

# Sertting psi=resb=0.0 because Sref doesnt require the adjoint solution
resb = np.zeros(n_vort)

llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
clb, cd, cdb, l, lb, d, db)

dSrefdX = xb
dSref0dalpha0 = alpha0b
dSrefdc = chordsb
rSref = gamab
print('dSrefdX', dSrefdX)
print('dSref0dalpha0',dSref0dalpha0)
print('dSrefdc', dSrefdc)
print('rSref', rSref)

###############################################################################
