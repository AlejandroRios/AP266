import numpy as np
from scipy.optimize import minimize

# Define Rosenbrock’s function
def objective(x):
 return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

exec(open('./rosen_grad.py').read())

# define number os design variables
dim_num = 2
dim_num = np.ones(dim_num)
# define initial point
x0 = dim_num*0
print(x0)
print(objective(x0))
nfev_BFGS = []

# Defining BFGS optimization problem
result =  minimize(objective, x0, method='BFGS', jac=rosen_grad,tol= 1e-14, options={'maxiter':100000, disp': True})


# Print results
print(result)
print(result.x)
print(result.fun)

nfev_BFGS.append(result.nfev + result.njev)
print(nfev_BFGS)
