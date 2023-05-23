from sympy import *


# these functions assume knowledge of all initial conditions
def ab2(y1_p, var1_name, var2_name, a, b, n, ics1, ics2):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = []
    for i in range(0, n + 1):
        ti.append(a + i * h)
    wi = [None] * (n+1)
    wi[0] = ics1
    wi[1] = ics2

    w_arr = [wi[0], wi[1]]
    for i in range(1, n):
        wi[i+1] = float(wi[i] + (h/2) * (3*y1_p.subs({v1: ti[i], v2: wi[i]}) - y1_p.subs({v1: ti[i-1], v2: wi[i-1]})))
        w_arr.append(wi[i+1])

    return ti, w_arr


def ab3(y1_p, var1_name, var2_name, a, b, n, ics1, ics2, ics3):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = []
    for i in range(0, n + 1):
        ti.append(a + i * h)
    wi = [None] * (n+1)
    wi[0] = ics1
    wi[1] = ics2
    wi[2] = ics3
    w_arr = [wi[0], wi[1], wi[2]]
    for i in range(2, n):
        wi[i+1] = float(wi[i] + (h/12) * (23*y1_p.subs({v1: ti[i], v2: wi[i]})
                                          - 16*y1_p.subs({v1: ti[i-1], v2: wi[i-1]})
                                          + 5*y1_p.subs({v1: ti[i-2], v2: wi[i-2]})))
        w_arr.append(wi[i+1])

    return ti, w_arr


def ab4(y1_p, var1_name, var2_name, a, b, n, ics1, ics2, ics3, ics4):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = []
    for i in range(0, n + 1):
        ti.append(a + i * h)
    wi = [None] * (n+1)
    wi[0] = ics1
    wi[1] = ics2
    wi[2] = ics3
    wi[3] = ics4
    w_arr = [wi[0], wi[1], wi[2], wi[3]]
    for i in range(3, n):
        wi[i+1] = float(wi[i] + (h/24) * (55*y1_p.subs({v1: ti[i], v2: wi[i]})
                                          - 59*y1_p.subs({v1: ti[i-1], v2: wi[i-1]})
                                          + 37*y1_p.subs({v1: ti[i-2], v2: wi[i-2]})
                                          - 9*y1_p.subs({v1: ti[i-3], v2: wi[i-3]})))
        w_arr.append(wi[i+1])

    return ti, w_arr


def ab5(y1_p, var1_name, var2_name, a, b, n, ics1, ics2, ics3, ics4, ics5):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = []
    for i in range(0, n + 1):
        ti.append(a + i * h)
    wi = [None] * (n+1)
    wi[0] = ics1
    wi[1] = ics2
    wi[2] = ics3
    wi[3] = ics4
    wi[4] = ics5
    w_arr = [wi[0], wi[1], wi[2], wi[3], wi[4]]
    for i in range(4, n):
        wi[i+1] = float(wi[i] + (h/720) * (1901*y1_p.subs({v1: ti[i], v2: wi[i]})
                                          - 2774*y1_p.subs({v1: ti[i-1], v2: wi[i-1]})
                                          + 2616*y1_p.subs({v1: ti[i-2], v2: wi[i-2]})
                                          - 1274*y1_p.subs({v1: ti[i-3], v2: wi[i-3]})
                                          + 251*y1_p.subs({v1: ti[i-4], v2: wi[i-4]})))
        w_arr.append(wi[i+1])

    return ti, w_arr


#################################################
var1 = 't'
var2 = 'y'
t = symbols(var1)
y = symbols(var2)
#################################################
y1_p_1 = t*E**(3*t) - 2 * y
sol = (1/5)*t*E**(3*t) - (1/25)*E**(3*t) + (1/25)*E**(-2*t)  # given solution
a1 = 0
b1 = 1
n1 = 5
ic1 = 0
ic2 = sol.subs(t, 0.2)
ic3 = sol.subs(t, 0.4)
ic4 = sol.subs(t, 0.6)
ic5 = sol.subs(t, 0.8)
t1, y1 = ab2(y1_p_1, var1, var2, a1, b1, n1, ic1, ic2)
print("Q1-a-AB2")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
t1, y1 = ab3(y1_p_1, var1, var2, a1, b1, n1, ic1, ic2, ic3)
print("Q1-a-AB3")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
t1, y1 = ab4(y1_p_1, var1, var2, a1, b1, n1, ic1, ic2, ic3, ic4)
print("Q1-a-AB4")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
t1, y1 = ab5(y1_p_1, var1, var2, a1, b1, n1, ic1, ic2, ic3, ic4, ic5)
print("Q1-a-AB5")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")