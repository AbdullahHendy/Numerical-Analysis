import numpy as np


def newton_divided_coeff(xi, yi):
    F = np.zeros((len(xi), len(yi)))  # create n x n empty array
    for i in range(0, len(xi)):  # initialize the first column of the table
        F[i][0] = yi[i]
    for i in range(1, len(xi)):
        for j in range(1, i + 1):
            F[i][j] = (F[i][j-1]-F[i-1][j-1])/(x[i]-x[i-j])

    coeff = []
    for k in range(0, len(xi)):
        coeff.append(F[k][k])
    return coeff


x = [0, 0.1, 0.3, 0.6, 1]
y = [-6, -5.89483, -5.65014, -5.17788, -4.28172]
ai = newton_divided_coeff(x, y)
print(ai)


