import numpy as np
from sympy import *


def neville_eval(xi, yi, x_appr):
    Q = np.zeros((len(xi), len(yi)))  # create n x n empty array
    for i in range(0, len(xi)):  # initialize the first column of the table
        Q[i][0] = yi[i]
    for i in range(1, len(xi)):
        for j in range(1, i + 1):
            Q[i][j] = ((x_appr - xi[i - j]) * Q[i, j - 1] - (x_appr - xi[i]) * Q[i - 1, j - 1]) / (xi[i] - xi[i - j])

    return Q[len(xi)-1, len(yi)-1]


x = [-0.75, -0.5, -0.25, 0]
y = [-0.0718125, -0.02475, 0.3349375, 1.101]
sol = neville_eval(x, y, 0.5)
print(sol)
