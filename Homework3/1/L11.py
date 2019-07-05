'''
Ex1 Homework3
function definition file
Author = Alejandro Rios
AeroStructural Optimization
ITA
'''
# Importing modules
import numpy as np
import matplotlib.pyplot as plt

# Define function to be tested
def f(x):
    fun = np.exp(x)/((np.cos(x))**3 + np.sin(x)**3)
    return fun
# function first derivative
def df(x):
    df_a = ((np.exp(x)*(np.cos(3*x) + np.sin(3*x)/2 + (3*np.sin(x))/2)) /
            (np.cos(x)**3 + np.sin(x)**3)**2)
    return df_a
# complex step derivation
def df_cs(x,h):
    df_cs = np.imag(f(x + 1j*h))/h
    return df_cs
# forward diff. derivation
def df_fd(x,h):
    df_fd = (f(x+h) - f(x))/h
    return df_fd

# point to test function
p_x = np.pi/4
# vector of h to be tested
h_v = np.logspace(-30,-1,20, endpoint=True)

error_cs = np.zeros(h_v.shape)
error_df = np.zeros(h_v.shape)

# Evaluating error with respect to exact derivative
for i, h in enumerate(h_v):
    error_cs[i] = np.abs(df_cs(p_x, h) - df(p_x))/np.abs(df(p_x))
    error_df[i] = np.abs(df_fd(p_x, h) - df(p_x))/np.abs(df(p_x))

# Generating Plot
plt.figure()
plt.loglog(h_v, error_cs,'-ob')
plt.loglog(h_v, error_df,'-ok')
plt.gca().invert_xaxis()
plt.legend(('CS', 'FD'),loc='center right')
plt.ylabel('Error')
plt.xlabel('h')
plt.grid(True)
plt.show()
