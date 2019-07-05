'''
Ex2.6 Homework3
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
from scipy.optimize import minimize
from scipy.optimize import root
import matplotlib.pyplot as plt

################################################################################
# EXECUTION

# Define number of vortex and other variables
n_vort = 200
b = 8.0
chord = 1.0
cla = 6.283

# Define kynks position
x = np.zeros([3,n_vort+1])
x[0] = x[0,:]
x[1] =  np.linspace(0, b, n_vort+1)
x[2] = x[2,:]
x = np.array(x,order='F')
# Define gama, alpha0, chords, cl0 local and cla local
gama = np.ones([n_vort])*2
gama[0] = 1
gama[-1] = 1
alpha0 = np.zeros([n_vort])
chords = np.ones([n_vort])*chord
cl0 = np.zeros([n_vort])
cla = np.ones([n_vort])*cla

# Define ângle of attack
alpha = 5.498*np.pi/180
# Define air desnity
rho = 1
# Define velocity component vector
vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
################################################################################
# Define input seed for reversal AD
resb = np.ones([n_vort])
srefb = -0.3
clb = 1.3
cdb = -0.5
lb = 3.0
db = 0.02

gama0 = np.zeros([n_vort])
sref, _, _, l, d = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)
###############################################################################
# Define Objective function
###############################################################################
def objfun(alpha0):
    # Define function that computes residuals for a given gama
    def resfunc(gama):
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res

    # Solve AED problem
    sol = root(resfunc,gama0)
    # Get solved gama values
    gama = sol.x
    # Compute CD function for the correct gama
    _, _, cd, _, _ = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)
    return cd

###############################################################################
# Define gradient of Objective function
###############################################################################

def objfungrad(alpha0):
    # Define function that compute residuals for given gama
    def resfunc(gama):
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res

    # Solve the AED problem
    sol = root(resfunc,gama0)
    # Get solved gama values
    gama = sol.x
    sref, _, _, l, d = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)
    # Initialize output variables for reverse calls
    cl = 0.0
    cd = 0.0
    res = np.zeros(n_vort)

    # Define function to compute the residual of the adjoint equation
    def adjfunbc(resb):
        # Define seeds
        srefb = 0.0
        clb = 0.0
        cdb = 1.0
        lb = 0.0
        db = 0.0
        xb = np.zeros([3,n_vort+1])
        xb = np.array(xb,order='F')
        gamab = np.zeros([n_vort])
        alpha0b = np.zeros([n_vort])
        chordsb = np.zeros([n_vort])

        # Run reverse mode
        llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
        chords, chordsb, cl0, cla, vinf, rho, res, resb.copy(), sref, srefb, cl,
        clb, cd, cdb, l, lb, d, db)
        # Compute the residual of the adjoint equation for CD
        adj_res = gamab
        return adj_res

    # Solve the adjoint problem
    psi0 = np.ones(n_vort)*10**6
    sol = root(adjfunbc, psi0)
    # Get adjoint variables
    psi = sol.x

    # Initialize derivative seeds
    srefb = 0
    clb = 0
    cdb = 1
    lb = 0
    db = 0
    xb = np.zeros([3,n_vort+1])
    xb = np.array(xb,order='F')

    gamab = np.zeros([1,n_vort])
    alpha0b = np.zeros([1,n_vort])
    chordsb = np.zeros([1,n_vort])
    resb = psi.copy()

    # Ruin the reverse mode
    llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
    chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
    clb, cd, cdb, l, lb, d, db)
    # Get the total ferivative
    dCDdalpha = alpha0b

    return dCDdalpha

###############################################################################
# Define constrain function
###############################################################################
def confun(alpha0):
    # Define function that computes residuals for a given gama
    def resfunc(gama) :
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res
    # Solve AED problem
    sol = root(resfunc,gama0)
    # Get solved gama values
    gama = sol.x
    # Compute CL function for the correct gama
    _, cl, _, _, _ = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)

    # Define equal constrain
    con = 0.5-cl
    return con

###############################################################################
# Define gradient of the constrain function
###############################################################################
def confungrad(alpha0):
    # Define function that compute residuals for given gama
    def resfunc(gama):
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res

    # Solve the AED problem
    sol = root(resfunc,gama0)
    # Get solved gama values
    gama = sol.x
    # Initialize output variables for reverse calls
    cl = 0.0
    cd = 0.0
    res = np.zeros(n_vort)

    # Define function to compute the residual of the adjoint equation
    def adjfunbc(resb):
        # Define seeds
        srefb = 0.0
        clb = 1.0
        cdb = 0.0
        lb = 0.0
        db = 0.0
        xb = np.zeros([3,n_vort+1])
        xb = np.array(xb,order='F')
        gamab = np.zeros([n_vort])
        alpha0b = np.zeros([n_vort])
        chordsb = np.zeros([n_vort])
        # Run reverse mode
        llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
        chords, chordsb, cl0, cla, vinf, rho, res, resb.copy(), sref, srefb, cl,
        clb, cd, cdb, l, lb, d, db)
        # Compute the residual of the adjoint equation for CD
        adj_res = gamab
        return adj_res

    # Solve the adjoint problem
    psi0 = np.ones(n_vort)*10**6
    sol = root(adjfunbc, psi0)

    # Get adjoint variables
    psi = sol.x
    # Initialize derivative seeds
    srefb = 0.0
    clb = 1.0
    cdb = 0.0
    lb = 0.0
    db = 0.0
    cl = 0.0
    cd = 0.0
    xb = np.zeros([3,n_vort+1])
    xb = np.array(xb,order='F')

    gamab = np.zeros([n_vort])
    alpha0b = np.zeros([n_vort])
    chordsb = np.zeros([n_vort])
    resb = psi.copy()

    # Run the reverse mode
    llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
    chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
    clb, cd, cdb, l, lb, d, db)

    # Get the total ferivative
    dCLdalpha0 = alpha0b

    congrad = -dCLdalpha0
    return congrad

################################################################################
# Otimization
# Define constrains type
con1 = {'type':'eq','fun':confun,'jac':confungrad}
# Run optimization
result = minimize (objfun, alpha0, jac = objfungrad, constraints = [con1], method = 'SLSQP',options={'maxiter': 100000, 'ftol': 1e-12, 'disp': False})

# Obtain optimal alpha0, CD and CL
alpha0 = result.x
CD = objfun(alpha0)
CL = 0.5-confun(alpha0)

def resfunc(gama):
    res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
    return res

sol = root(resfunc,gama0)
gama = sol.x

print(CL)
print(CD)
print(gama)

print('CL:',CL)
print('CD:',CD)
print('gama:',gama)
print('alpha0:', alpha0*180/np.pi)

y1 = x[1,:-1]
y2 = x[1,1:]
y = (y2+y1)/2


##########################################################################################

vinf = vinf[0]
CL_ell = 0.5

#yelli = np.linspace(-b/2,b/2,16)
yelli = np.linspace(0,b,n_vort)
Gama0 = ((2.0*vinf*sref)/(np.pi*b))*CL
Gama_ell = 2.0*Gama0 *np.sqrt((yelli/b)*(1.0-(yelli/b)))
CD_ell= (CL**2 * sref)/(np.pi*b**2)

print(CD)
#Gama_ell = Gama0 *np.sqrt(1.0-(2*yelli/b)**2)


##########################################################################################
plt.figure(1)
plt.plot(y,gama,'--b',yelli,Gama_ell,'--k')
plt.title("")
plt.xlabel("b [units]")
plt.ylabel("gama")
# plt.grid(True)
plt.legend(('Optimized', 'Elliptical'),loc='north right')

plt.figure(2)
plt.plot(y,alpha0*180/np.pi,'-b')
plt.title("")
plt.xlabel("b [units]")
plt.ylabel("alpha0 [deg]")
# plt.grid(True)

plt.figure(3)
plt.semilogy(n_vort,CL,'--ob',n_vort,CL_ell,'--ob')
plt.title("")
plt.xlabel("n_vort")
plt.ylabel("CL")
plt.grid(True)

plt.figure(4)
plt.plot(n_vort,CD,'--ob',n_vort,CD_ell,'--ok')
plt.title("")
plt.xlabel("n_vort")
plt.ylabel("CD")
plt.grid(True)

plt.show()
