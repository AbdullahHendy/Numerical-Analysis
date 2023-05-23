from sympy import *


def gen_n_aitken(func, var_name, p0, iter_lim):
    v = symbols(var_name)
    p0 = p0
    p1 = func.subs(v, p0)
    p2 = func.subs(v, p1)
    for i in range(0, iter_lim + 1):
        p_hat = float(p0 - (((p1 - p0) ** 2) / (p2 - 2 * p1 + p0)))
        print(f"P{i} = {p_hat}")
        p0 = p1
        p1 = p2
        p2 = func.subs(v, p2)


def newton(func, var_name, p_0, tolerance, iter_lim):
    v = symbols(var_name)
    func_prime = diff(func)
    for i in range(1, iter_lim + 1):
        p = float(p_0 - func.subs(v, p_0) / func_prime.subs(v, p_0))
        if p == p_0:
            print(f"{p} is an exact root found after {i} iterations")
            return
        elif abs(p - p_0) < tolerance:
            print(f"{p} is a root found after {i} iterations with absolute error on the order of {tolerance}")
            return
        p_0 = p
    print(f"failed to find a root in {iter_lim} iterations")


def newton_aitken(func, var_name, p0, tolerance, iter_lim):
    v = symbols(var_name)
    p0 = p0
    p1 = float(p0 - func.subs(v, p0) / diff(func).subs(v, p0))
    p2 = float(p1 - func.subs(v, p1) / diff(func).subs(v, p1))
    p_hat_1 = 2 ** 30  # place holder

    for i in range(0, iter_lim + 1):
        p_hat = float(p0 - (((p1 - p0) ** 2) / (p2 - 2 * p1 + p0)))
        p0 = p1
        p1 = p2
        p2 = float(p2 - func.subs(v, p2) / diff(func).subs(v, p2))
        if abs(p_hat - p_hat_1) < tolerance:
            print(f"To get absolute error within {tolerance}, Aitken's method gives p{i} = {p_hat}")
            return
        p_hat_1 = p_hat
    print("Aitken's method failed")


def q13_aitken(func, var_name, iter_lim):
    v = symbols(var_name)
    p1 = func.subs(v, 1)
    p2 = func.subs(v, 2)
    p3 = func.subs(v, 3)
    j = 4
    for i in range(1, iter_lim + 1):
        p_hat = float(p1 - (((p2 - p1) ** 2) / (p3 - 2 * p2 + p1)))
        p1 = p2
        p2 = p3
        p3 = func.subs(v, j)
        if abs(p_hat) < 5*10**(-2):
            print(f"Aitken's method generates p{i} = {p_hat} that is less than {5*10**(-2)} margin")
            return
        j = j + 1
    print("Failed")


var = 'x'
x = symbols(var)
tol = 5 * 10 ** (-2)
# expression = (2 - E ** x + x ** 2) / 3
# expression = ((E ** x) / 3) ** (1 / 2)
# expression = E ** (6 * x) + 3 * (ln(2) ** 2) * E ** (2 * x) - ln(8) * E ** (4 * x) - (ln(2) ** 3)
p_init_0 = 0.75
N = 100

# gen_n_aitken(expression, var, p_init_0, N)
# newton(expression, var, p_init_0, tol, N)
# newton_aitken(expression, var, p_init_0, tol, N)
#
q13_a = 1 / x
q13_b = 1 / (x**2)
q13_aitken(q13_b, var, N)
