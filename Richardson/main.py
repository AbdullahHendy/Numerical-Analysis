import numpy as np
from sympy import *


# This function only works when M = N_1(h) + k1 h^2 + k2 h^4 + k3 h^6 + ....
def richardson_extrapolation_given(ni, n):  # n represents N_n(h)
    f = np.zeros((n, n))  # create n x n empty array
    for i in range(0, n):  # initialize the first column of the table given approximations (data)
        f[i][0] = ni[i]
    for i in range(1, n):
        for j in range(1, i + 1):
            f[i][j] = ((4 ** j) * f[i][j - 1] - f[i - 1][j - 1]) / ((4 ** j) - 1)

    return f[n - 1][n - 1]


# This function only works when M = N_1(h) + k1 h^2 + k2 h^4 + k3 h^6 + ....
def richardson_extrapolation_central_diff(func, var_name, h, n, x_appr):  # n represents N_n(h)
    f = np.zeros((n, n))  # create n x n empty array
    v = symbols(var_name)
    for i in range(0, n):  # initialize the first column of the table given approximations (data)
        f[i][0] = (func.subs(v, x_appr + (h / (2 ** i))) - func.subs(v, x_appr - (h / (2 ** i)))) / (2 * (h / (2 ** i)))
    for i in range(1, n):
        for j in range(1, i + 1):
            f[i][j] = ((4 ** j) * f[i][j - 1] - f[i - 1][j - 1]) / ((4 ** j) - 1)

    return f[n - 1][n - 1]


#################################################
n_i = [1.570796, 1.896119, 1.974232, 1.99357]
order_1 = 4
n_4 = richardson_extrapolation_given(n_i, order_1)
print(n_4)
#################################################

var = 'x'
x = symbols(var)
expression = x * (E ** x)
h_2 = 0.2
order_2 = 3
x_approximate = 2

n_3 = richardson_extrapolation_central_diff(expression, var, h_2, order_2, x_approximate)
print(n_3)
