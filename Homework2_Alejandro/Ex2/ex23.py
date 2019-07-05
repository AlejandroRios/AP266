import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution

# Define Rosenbrock’s function
def objective(x):
 return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

exec(open('./rosen_grad.py').read())




# Initialize list of number of function eval.
nfev_NM = []
ofv_NM = []
nfev_DE = []
ofv_DE = []
nfev_BFGS = []
ofv_BFGS = []
nfev_CG = []
ofv_CG = []
nv = []
bounds=[]
for i in range(2,21):
    # define number os design variables
    dim_num = i
    dim_num = np.ones(dim_num)
    # define initial point
    x0 = dim_num*0

    seed=np.random.seed(seed=3)
    # define number os design variables
    dim_num = i
    dim_num = np.ones(dim_num)
# define initial point
    x0 = dim_num*0
# Defining bounds
    bv = 4
    bounds.append((-bv,bv))
    # Defining bounds
    bv = 2
    bounds.append((-bv,bv))
    # Defining Genetic Algorithm problem
    # result0 = differential_evolution(objective, bounds, maxiter=100000, polish=False, options={'disp': True})
    result0 = differential_evolution(objective, bounds,maxiter=100000,polish=False)
    # Defining Nelder-Mead optimization problem
    result1 = minimize(objective, x0, method='Nelder-Mead', tol=1e-14)

    result2 =   minimize(objective, x0, method='CG', jac=rosen_grad,tol=1e-14, options={'maxiter':10000})

    result3 =  minimize(objective, x0, method='BFGS', jac=rosen_grad,tol= 1e-14, options={'maxiter':10000})
    nfev_DE.append(result0.nfev)
    nfev_NM.append(result1.nfev)
    nfev_BFGS.append(result2.nfev + result2.njev)
    nfev_CG.append(result3.nfev + result3.njev)
    nv.append(i)

    ofv_DE.append(result0.fun)
    ofv_NM.append(result1.fun)
    ofv_BFGS.append(result2.fun)
    ofv_CG.append(result3.fun)
    # nv.append(i)

    plt.figure(1)
    # plt.semilogy(nv,nfev_NM,'-ob',nv,nfev_DE,'-om',nv,nfev_BFGS,'-or',nv,nfev_CG,'-oy')
    plt.plot(nv,nfev_NM,'-ob',nv,nfev_DE,'-om',nv,nfev_BFGS,'-or',nv,nfev_CG,'-oy')
    plt.figure(2)
    # plt.semilogy(nv,ofv_NM,'-ob',nv,ofv_DE,'-om',nv,ofv_BFGS,'-or',nv,ofv_CG,'-oy')
    plt.plot(nv,ofv_NM,'-ob',nv,ofv_DE,'-om',nv,ofv_BFGS,'-or',nv,ofv_CG,'-oy')


plt.figure(1)
plt.title("NM, BFGS and CG method comparision")
plt.xlabel("n. design variables")
plt.ylabel("Num. fun. eval. ")
plt.grid(True)
plt.legend(('NM','DE', 'BFGS', 'CG'),
           loc='upper left')

plt.figure(2)
plt.title("NM,DE, BFGS and CG method comparision")
plt.xlabel("n. design variables")
plt.ylabel("Obj. fun. value ")
plt.grid(True)
plt.legend(('NM','DE','BFGS', 'CG'),
           loc='upper left')

plt.show()
