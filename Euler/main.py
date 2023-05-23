from sympy import *


def euler_method(func, var1_name, var2_name, a, b, n, ic):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = a
    wi = ic
    t_arr = [ti]
    w_arr = [wi]
    for i in range(1, n + 1):
        wi = float(wi + h * func.subs({v1: ti, v2: wi}))
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
expression1 = 1 + y / t
a1 = 1
b1 = 2
n1 = 4
ic1 = 2
t1, y1 = euler_method(expression1, var1, var2, a1, b1, n1, ic1)
print("Q1-c")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
#################################################
expression2 = cos(2 * t) + sin(3 * t)
a2 = 0
b2 = 1
n2 = 4
ic2 = 1
t2, y2 = euler_method(expression2, var1, var2, a2, b2, n2, ic2)
print("Q1-d")
[print(i, end=' ') for i in t2]
print("")
[print(i, end=' ') for i in y2]
print("\n---------------------------------------------------")
#################################################
expression3 = 2*(y/t) + (t**2) * (E**t)
a3 = 1
b3 = 2
n3 = 10
ic3 = 0
t3, y3 = euler_method(expression3, var1, var2, a3, b3, n3, ic3)
print("Q9-a")
[print(i, end=' ') for i in t3]
print("")
[print(i, end=' ') for i in y3]
print("\n---------------------------------------------------")
#################################################
expression4 = -1*y + t + 1
a4 = 0
b4 = 5
n4_1 = 25
ic4 = 1
t4_1, y4_1 = euler_method(expression4, var1, var2, a4, b4, n4_1, ic4)
print("Q11-a - h=0.2")
print(y4_1[n4_1])
print("Q11-a - h=0.1")
n4_2 = 50
t4_2, y4_2 = euler_method(expression4, var1, var2, a4, b4, n4_2, ic4)
print(y4_2[n4_2])
print("Q11-a - h=0.05")
n4_3 = 100
t4_3, y4_3 = euler_method(expression4, var1, var2, a4, b4, n4_3, ic4)
print(y4_3[n4_3])
print("\n---------------------------------------------------")
#################################################
