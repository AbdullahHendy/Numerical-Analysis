import numpy as np
from sympy import *


def hermite_coeff(xi, yi, ypi):
    n = len(xi) - 1
    z = np.zeros(2 * n + 2)
    q = np.zeros((2 * n + 2, 2 * n + 2))  # create n x n empty array

    for i in range(0, n + 1):
        z[2 * i] = xi[i]
        z[2 * i + 1] = xi[i]
        q[2 * i][0] = yi[i]
        q[2 * i + 1][0] = yi[i]
        q[2 * i + 1][1] = ypi[i]
        if i != 0:
            q[2 * i][1] = (q[2 * i][0] - q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1])
    for i in range(2, 2 * n + 2):
        for j in range(2, i + 1):
            q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / (z[i] - z[i - j])

    coeff = []
    for k in range(0, 2 * n + 2):
        coeff.append(q[k][k])
    return coeff


def hermite_3(xi, coeff, x_appr):  # hermite polynomial of degree 3
    h_3 = 0
    vari = 'x'
    x = symbols(vari)
    h = 1
    n = 2
    for i in range(0, 2 * n):
        if i < 3:
            h = (x - xi[0]) ** i
        elif i < 5:
            h = h * (x - xi[1])

        h_3 = h_3 + h * coeff[i]
    return h_3.subs(x, x_appr)


def hermite_5(xi, coeff, x_appr):  # hermite polynomial of degree 5
    h_5 = 0
    vari = 'x'
    x = symbols(vari)
    h = 1
    n = 3
    for i in range(0, 2 * n):
        if i < 3:
            h = (x - xi[0]) ** i
        elif i < 5:
            h = h * (x - xi[1])
        elif i < 7:
            h = h * (x - xi[2])

        h_5 = h_5 + h * coeff[i]
    return h_5.subs(x, x_appr)


def hermite_7(xi, coeff, x_appr):  # hermite polynomial of degree 7
    h_7 = 0
    vari = 'x'
    x = symbols(vari)
    h = 1
    n = 4
    for i in range(0, 2 * n):
        if i < 3:
            h = (x - xi[0]) ** i
        elif i < 5:
            h = h * (x - xi[1])
        elif i < 7:
            h = h * (x - xi[2])
        elif i < 9:
            h = h * (x - xi[3])
        h_7 = h_7 + h * coeff[i]
    return h_7.subs(x, x_appr)


####################################################################
x_1a = [8.3, 8.6]
y_1a = [17.56492, 18.50515]
yp_1a = [3.116256, 3.151762]
x_appr_1a = 8.4
coe_1a = hermite_coeff(x_1a, y_1a, yp_1a)
print(f"hermite coefficient: {coe_1a} ")
value_1a = hermite_3(x_1a, coe_1a, x_appr_1a)
print(f"the approximation of {x_appr_1a} using 3rd degree hermite polynomial is {value_1a}")
print("------------------------------------------------------------------------------------------")
####################################################################
x_1c = [-0.50, -0.25, 0]
y_1c = [-0.02475, 0.3349375, 1.101]
yp_1c = [0.751, 2.189, 4.002]
x_appr_1c = -(1 / 3)
coe_1c = hermite_coeff(x_1c, y_1c, yp_1c)
print(f"hermite coefficient: {coe_1c} ")
value_1c = hermite_3(x_1c, coe_1c, x_appr_1c)
print(f"the approximation of {x_appr_1c} using 3rd degree hermite polynomial is {value_1c}")
print("------------------------------------------------------------------------------------------")
####################################################################
x_5a = [0.3, 0.32, 0.35]
y_5a = [0.29552, 0.31457, 0.3429]
yp_5a = [0.95534, 0.94924, 0.93937]
x_appr_5a = 0.34
coe_5a = hermite_coeff(x_5a, y_5a, yp_5a)
print(f"hermite coefficient: {coe_5a} ")
value_5a = hermite_5(x_5a, coe_5a, x_appr_5a)
print(f"the approximation of {x_appr_5a} using 5th degree hermite polynomial is {value_5a}")
print("------------------------------------------------------------------------------------------")
####################################################################
x_5c = [0.3, 0.32, 0.33, 0.35]
y_5c = [0.29552, 0.31457, 0.32404, 0.3429]
yp_5c = [0.95534, 0.94924, 0.94604, 0.93937]
x_appr_5c = 0.34
coe_5c = hermite_coeff(x_5c, y_5c, yp_5c)
print(f"hermite coefficient: {coe_5c} ")
value_5c = hermite_7(x_5c, coe_5c, x_appr_5c)
print(f"the approximation of {x_appr_5c} using 7th degree hermite polynomial is {value_5c}")
