from sympy import *


# n must be even
def composite_simpson(func, var_name, a, b, n):
    v = symbols(var_name)
    h = (b - a) / n
    first_term = func.subs(v, a)
    last_term = func.subs(v, b)
    sum_odd = 0
    sum_even = 0
    for i in range(1, n):
        s = a + i * h
        if i % 2 == 0:
            sum_even = sum_even + func.subs(v, s)
        else:
            sum_odd = sum_odd + func.subs(v, s)

    final_approx = float((h / 3) * (first_term + 2 * sum_even + 4 * sum_odd + last_term))
    return final_approx


# n must be even
def composite_trapezoidal(func, var_name, a, b, n):
    v = symbols(var_name)
    h = (b - a) / n
    first_term = func.subs(v, a)
    last_term = func.subs(v, b)
    sum_all = 0
    for i in range(1, n):
        s = a + i * h
        sum_all = sum_all + func.subs(v, s)

    final_approx = float((h / 2) * (first_term + 2 * sum_all + last_term))
    return final_approx


#################################################
var = 'x'
x = symbols(var)
#################################################
expression1 = x * ln(x)
a1 = 1
b1 = 2
n1 = 4
integral1_simp = composite_simpson(expression1, var, a1, b1, n1)
integral1_trap = composite_trapezoidal(expression1, var, a1, b1, n1)
print(f"Q3-a: {integral1_simp}")
print(f"Q1-a: {integral1_trap}")
#################################################
expression2 = x ** 2 * cos(x)
a2 = 0
b2 = pi
n2 = 6
integral2_simp = composite_simpson(expression2, var, a2, b2, n2)
integral2_trap = composite_trapezoidal(expression2, var, a2, b2, n2)
print(f"Q3-d: {integral2_simp}")
print(f"Q1-d: {integral2_trap}")
#################################################
expression3 = (x ** 2) * E**(-1*x**2)
a3 = 0
b3 = 2
n3 = 8
integral3_trap = composite_trapezoidal(expression3, var, a3, b3, n3)
integral3_simp = composite_simpson(expression3, var, a3, b3, n3)
print(f"Q8-a: {integral3_trap}")
print(f"Q8-b: {integral3_simp}")
#################################################
