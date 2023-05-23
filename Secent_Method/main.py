from sympy import *


def secent(func, var_name, p_0, p_1, tolerance, iter_lim):
    v = symbols(var_name)
    fp_0 = func.subs(v, p_0)
    fp_1 = func.subs(v, p_1)
    for i in range(1, iter_lim + 1):
        p = float(p_1 - (fp_1 * (p_1 - p_0)) / (fp_1 - fp_0))
        if p == p_1:
            print(f"exact root found after {i} iterations")
            return p
        elif abs(p - p_1) / abs(p) < tolerance:
            print(f"a root found after {i} iterations with relative error on the order of {tolerance}")
            return p
        p_0 = p_1
        p_1 = p
        fp_0 = fp_1
        fp_1 = func.subs(v, p)
    return f"failed to find a root in {iter_lim} iterations"


var = 'x'
x = symbols(var)
expression = (x-2)**2 - ln(x)
tol = 10 ** (-5)
p_init_0 = E
p_init_1 = 4
N = 1000

secent_p = secent(expression, var, p_init_0, p_init_1, tol, N)
print(secent_p)
