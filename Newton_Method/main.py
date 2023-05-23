from sympy import *


def newton(func, var_name, p_0, tolerance, iter_lim):
    v = symbols(var_name)
    func_prime = diff(func)
    for i in range(1, iter_lim + 1):
        p = float(p_0 - func.subs(v, p_0) / func_prime.subs(v, p_0))
        if p == p_0:
            print(f"exact root found after {i} iterations")
            return p
        elif abs(p - p_0) / abs(p) < tolerance: # relative error is more conservative than absolute error
            print(f"a root found after {i} iterations with relative error on the order of {tolerance}")
            return p
        p_0 = p
    return f"failed to find a root in {iter_lim} iterations"


var = 'x'
x = symbols(var)
expression = E ** (6 * x) + 3 * (ln(2) ** 2) * E ** (2 * x) - ln(8) * E ** (4 * x) - (ln(2) ** 3)
tol = 2 * 10 ** (-4)
p_init = 0
N = 1000

newton_p = newton(expression, var, p_init, tol, N)
print(newton_p)
