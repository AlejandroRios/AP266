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
n_vort = 16
b = 8
chord = 0.5
cla = 6.283



x = np.zeros([3,n_vort+1])
x[0] = x[0,:]
x[1] =  np.linspace(0, b, n_vort+1)
x[2] = x[2,:]
x = np.array(x,order='F')

gama = np.ones([n_vort])*2
gama[0] = 1
gama[-1] = 1


alpha0 = np.zeros([n_vort])

chords = np.ones([n_vort])*chord

cl0 = np.zeros([n_vort])

cla = np.ones([n_vort])*cla



alpha = 5.498*np.pi/180
rho = 1
vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])


resb = np.ones([n_vort])
srefb = -0.3
clb = 1.3
cdb = -0.5
lb = 3.0
db = 0.02

###############################################################################
# Define Objective function

def objfun(chords):
    gama0 = np.zeros([1,n_vort])
    def resfunc(gama):
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res

    sol = root(resfunc,gama0)
    gama = sol.x
    _, _, CD, _, _ = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)

    return CD

###############################################################################
# Define gradient of Objective function

def objfungrad(chords):
    gama0 = np.zeros([n_vort])
    def resfunc(gama):
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res

    sol = root(resfunc,gama0)
    gama = sol.x
    res, sref, cl, cd, l, d = llt.tapenade_main(x, gama, alpha0, chords, cl0, cla, vinf, rho)

    def adjfunbc(resb):

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

        llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
        chords, chordsb, cl0, cla, vinf, rho, res, resb.copy(), sref, srefb, cl,
        clb, cd, cdb, l, lb, d, db)

        adj_res = gamab
        return adj_res

    psi0 = np.ones(n_vort)*10**6
    sol = root(adjfunbc, psi0)

    # print(sol)

    psi = sol.x

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

    llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
    chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
    clb, cd, cdb, l, lb, d, db)

    dCDdchord = chordsb

    return dCDdchord

###############################################################################
# Define constrain function

def confun(chords):
    gama0 = np.zeros([n_vort])
    def resfunc(gama) :
     res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
     return res

    sol = root(resfunc,gama0)

    gama = sol.x
    res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)

    _, CL, _, _, _ = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)

    CLstr = 0.5
    con = CLstr -CL

    return con

# ###############################################################################
# Define gradient of the constrain function

def confungrad(chords):
    gama0 = np.zeros([n_vort])
    def resfunc(gama):
        res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
        return res

    sol = root(resfunc,gama0)
    gama = sol.x
    res, sref, cl, cd, l, d = llt.tapenade_main(x, gama, alpha0, chords, cl0, cla, vinf, rho)
    def adjfunbc(resb):

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

        llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
        chords, chordsb, cl0, cla, vinf, rho, res, resb.copy(), sref, srefb, cl,
        clb, cd, cdb, l, lb, d, db)

        adj_res = gamab
        return adj_res

    psi0 = np.ones(n_vort)*10**6
    sol = root(adjfunbc, psi0)

    # print(sol)

    psi = sol.x

    srefb = 0
    clb = 1
    cdb = 0
    lb = 0
    db = 0


    xb = np.zeros([3,n_vort+1])
    xb = np.array(xb,order='F')

    gamab = np.zeros([n_vort])
    alpha0b = np.zeros([n_vort])
    chordsb = np.zeros([n_vort])
    resb = psi.copy()

    llt_b.tapenade_main_b(x, xb, gama, gamab, alpha0, alpha0b,
    chords, chordsb, cl0, cla, vinf, rho, res, resb, sref, srefb, cl,
    clb, cd, cdb, l, lb, d, db)

    dCLdchord = chordsb
    congrad = -dCLdchord
    return congrad

################################################################################
# Otimization

con1 = {'type':'eq',
'fun':confun,
'jac':confungrad}

# bounds = [[0.0 , None ]]* n_vort


result = minimize (objfun, chords, jac = objfungrad, constraints = [con1], method = 'SLSQP',options={'maxiter': 100000, 'ftol': 1e-12,'disp': True})
gama0 = np.zeros([n_vort])
def resfunc(gama):
    res = llt.get_residuals(x, gama, alpha0, chords, cl0, cla, vinf, rho)
    return res

sol = root(resfunc,gama0)
gama = sol.x
chords = result.x

Sref, CL, CD, L, D = llt.get_functions(x, gama, alpha0, chords, cl0, cla, vinf, rho)

print('CL:',CL)
print('CD:',CD)
print('gama:',gama)
print('alpha0:', alpha0)
print('alpha0:', chords)

y1 = x[1,:-1]
y2 = x[1,1:]
y = (y2+y1)/2

#############################################################################
plt.figure(1)
plt.plot(y,gama,'--b')
plt.title("")
plt.xlabel("b")
plt.ylabel("gama")
plt.grid(True)
# plt.legend(('NM','DE','BFGS', 'CG'),
#            loc='upper left')

plt.figure(2)
plt.plot(y,chords,'--b')
plt.title("")
plt.xlabel("b")
plt.ylabel("alpha0")
plt.grid(True)

plt.show()
