'''
01_gd2d.py
uses the gradient descent to minimize a 2D function
'''

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import auxmod as am
from scipy.optimize import minimize_scalar
from scipy.optimize import line_search

# USER INPUTS

# Initialize function counter (we use a list so that we can change its value
# within objfun)
nfunc_calls = [0]

# Get the objective function
exec(open('./radfun.py').read())

# Number of iterations
maxiter = 100

# Line search flag
ls_option = None

# EXECUTION

# Create figure and axes handle
fig = plt.figure()
ax = plt.gca()

# Plot contour
am.plot_contour(objfun, ax,
                xmin, xmax, ymin, ymax, zmin, zmax)

# Initialize list of evaluated points and function values
xk = []
fk = []

# Set starting point
xk.append(x0)

# Compute optimum value
fopt = objfun(xopt)

# Reinitialize function counter since it was already used during the contour
# plots
nfunc_calls = [0]

# Iterate
for k in range(maxiter):

    # Evaluate objective function
    f = objfun(xk[-1])

    # Termination criteria
    if np.abs(f-fopt) < 1e-6:
        break

    # Evaluate gradient
    gradf = objfungrad(xk[-1])

    # Evaluate the hessian
    Hf = objfunhess(xk[-1])

    # Define search direction
    pk = np.linalg.solve(Hf,-gradf)

    ## LINE SEARCH
    if ls_option == 'exact':

        # Define univariate Function
        def unifun(alpha):
            xx = xk[-1] + alpha*pk
            ff = objfun(xx)
            return ff

        # Run optimizer to find alpha
        result = minimize_scalar(unifun,tol=1e-6)
        alpha = result.x

    elif ls_option == 'inexact':
        alpha,_,_,_,_,_ = line_search(objfun, objfungrad, xk[-1], pk,gfk=gradf, c2=0.1)

    elif ls_option == None:

        alpha = 1.0

    # Compute next iteration point
    xnext = xk[-1] + alpha*pk

    # Store results to the lists
    fk.append(f)
    xk.append(xnext)

# Print function calls
print('Iterations',k)
print('Function calls',nfunc_calls)

# Plot optimization history
am.plot_path(ax, xk, xopt)

# Plot the objective function history
am.plot_history(ax, fk, xopt, objfun)

plt.show()
