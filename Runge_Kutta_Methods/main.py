from sympy import *


def modified_euler(y1_p, var1_name, var2_name, a, b, n, ic):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = a
    wi = ic
    t_arr = [ti]
    w_arr = [wi]
    for i in range(1, n + 1):
        ti_1 = a + i * h  # ti+1
        wi = float(wi + (h/2) * (y1_p.subs({v1: ti, v2: wi})
                                 + y1_p.subs({v1: ti_1, v2: wi + h*y1_p.subs({v1: ti, v2: wi})})))
        ti = a + i * h
        t_arr.append(ti)
        w_arr.append(wi)

    return t_arr, w_arr


def mid_point(y1_p, var1_name, var2_name, a, b, n, ic):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = a
    wi = ic
    t_arr = [ti]
    w_arr = [wi]
    for i in range(1, n + 1):
        wi = float(wi + h * (y1_p.subs({v1: ti + (h/2), v2: wi + (h/2)*y1_p.subs({v1: ti, v2: wi})})))
        ti = a + i * h
        t_arr.append(ti)
        w_arr.append(wi)

    return t_arr, w_arr


def r_k_4(y1_p, var1_name, var2_name, a, b, n, ic):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = a
    wi = ic
    t_arr = [ti]
    w_arr = [wi]
    for i in range(1, n + 1):
        k1 = h * y1_p.subs({v1: ti, v2: wi})
        k2 = h * y1_p.subs({v1: ti + (h/2), v2: wi + (k1/2)})
        k3 = h * y1_p.subs({v1: ti + (h/2), v2: wi + (k2/2)})
        k4 = h * y1_p.subs({v1: ti + h, v2: wi + k3})

        wi = float(wi + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
        ti = a + i * h
        t_arr.append(ti)
        w_arr.append(wi)

    return t_arr, w_arr


#################################################
var1 = 't'
var2 = 'y'
t = symbols(var1)
y = symbols(var2)
#################################################
y1_prime_1 = 1 + (y/t) + (y/t)**2
a1 = 1
b1 = 3
n1 = 10
ic1 = 0
t1, y1 = modified_euler(y1_prime_1, var1, var2, a1, b1, n1, ic1)
print("Q3-b - modified euler")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
t1, y1 = mid_point(y1_prime_1, var1, var2, a1, b1, n1, ic1)
print("Q7-b - midpoint")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
t1, y1 = r_k_4(y1_prime_1, var1, var2, a1, b1, n1, ic1)
print("Q15-b - RK4")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
print("\n---------------------------------------------------")
#################################################
y1_prime_2 = -5*y + 5*t**2 + 2*t
a2 = 0
b2 = 1
n2 = 10
ic2 = 1/3
t2, y2 = modified_euler(y1_prime_2, var1, var2, a2, b2, n2, ic2)
print("Q3-d - modified euler")
[print(i, end=' ') for i in t2]
print("")
[print(i, end=' ') for i in y2]
print("\n---------------------------------------------------")
t2, y2 = mid_point(y1_prime_2, var1, var2, a2, b2, n2, ic2)
print("Q7-d - midpoint")
[print(i, end=' ') for i in t2]
print("")
[print(i, end=' ') for i in y2]
print("\n---------------------------------------------------")
t2, y2 = r_k_4(y1_prime_2, var1, var2, a2, b2, n2, ic2)
print("Q15-d - RK4")
[print(i, end=' ') for i in t2]
print("")
[print(i, end=' ') for i in y2]
print("\n---------------------------------------------------")
#################################################
