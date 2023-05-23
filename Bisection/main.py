from sympy import *


def bisection(func, var_name, lim_a, lim_b, tolerance, iter_lim):
    if abs(lim_a) == abs(lim_b):
        lim_a = lim_a - 1
    v = symbols(var_name)
    fa = func.subs(v, lim_a)
    # random large number to be updated after the first iteration.
    # sys.maxsize wasn't used to prevent possible overflow when checking for tolerance.
    p = 2 ** 30
    for i in range(1, iter_lim + 1):
        p_1 = p
        p = lim_a + (lim_b - lim_a) / 2
        fp = func.subs(v, p)
        if fp == 0:
            print(f"exact root found after {i} iterations")
            return p
        elif abs(p - p_1) / abs(p) < tolerance:
            print(f"a root found after {i} iterations with relative error on the order of {tolerance}")
            return p
        if fa * fp > 0:
            lim_a = p
            fa = fp
        else:
            lim_b = p
    return f"failed to find a root in {iter_lim} iterations"


var = 'x'
x = symbols(var)
expression = x - 2 * sin(x)
tol = 10 ** (-5)
a = 1.5
b = 3
N = 10000

root = bisection(expression, var, a, b, tol, N)
print(root)
