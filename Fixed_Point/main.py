from sympy import *


def fixed_point(func, var_name, p_0, tolerance, iter_lim):
    v = symbols(var_name)
    for i in range(1, iter_lim + 1):
        p = func.subs(v, p_0)
        if p == p_0:
            print(f"exact fixed point found after {i} iterations")
            return p
        elif abs(p - p_0) < tolerance:
            print(f"a fixed point found after {i} iterations with relative error on the order of {tolerance}")
            return p
        p_0 = p
    return f"failed to find a fixed point in {iter_lim} iterations"


var = 'x'
x = symbols(var)
expression = ((E**x)/3)**(1/2)
tol = 10 ** (-5)
p_init = 0.5
N = 10000

fixed_p = fixed_point(expression, var, p_init, tol, N)
print(fixed_p)
