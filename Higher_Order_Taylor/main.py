from sympy import *


# It would be better to only input the y prime of t and let the code do the high order derivatives but
# this implementation is chosen for the simplicity of the code
def taylor_two(y1_p, y2_p, var1_name, var2_name, a, b, n, ic):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = a
    wi = ic
    t_arr = [ti]
    w_arr = [wi]
    for i in range(1, n + 1):
        wi = float(wi + h * y1_p.subs({v1: ti, v2: wi}) + ((h ** 2) / 2) * y2_p.subs({v1: ti, v2: wi}))
        ti = a + i * h
        t_arr.append(ti)
        w_arr.append(wi)

    return t_arr, w_arr


def taylor_four(y1_p, y2_p, y3_p, y4_p, var1_name, var2_name, a, b, n, ic):
    v1 = symbols(var1_name)
    v2 = symbols(var2_name)
    h = (b - a) / n
    ti = a
    wi = ic
    t_arr = [ti]
    w_arr = [wi]
    for i in range(1, n + 1):
        wi = float(wi + h * y1_p.subs({v1: ti, v2: wi}) + ((h ** 2) / 2) * y2_p.subs({v1: ti, v2: wi})
                   + ((h ** 3) / 6) * y3_p.subs({v1: ti, v2: wi}) + ((h ** 4) / 24) * y4_p.subs({v1: ti, v2: wi}))
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
y_1prime_1 = E ** (t - y)
y_2prime_1 = (1 - y_1prime_1) * E ** (t - y)
y_3prime_1 = E ** (t - y) * ((1-y_1prime_1)**2 - y_2prime_1)
y_4prime_1 = E ** (t - y) * ((1-y_1prime_1)**3 - 3*y_2prime_1 + 3*y_1prime_1*y_2prime_1 - y_3prime_1)
a1 = 0
b1 = 1
n1 = 2
ic1 = 1
t1, y1 = taylor_two(y_1prime_1, y_2prime_1, var1, var2, a1, b1, n1, ic1)
print("Q2-a")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
t1, y1 = taylor_four(y_1prime_1, y_2prime_1, y_3prime_1, y_4prime_1, var1, var2, a1, b1, n1, ic1)
print("Q4-a")
[print(i, end=' ') for i in t1]
print("")
[print(i, end=' ') for i in y1]
print("\n---------------------------------------------------")
#################################################
y_1prime_2 = (1 + t) / (1 + y)
y_2prime_2 = ((1 + y) - (1 + t) * y_1prime_2) / (1 + y) ** 2
y_3prime_2 = (-2*y_1prime_2*(y+1-(t+1)*y_1prime_2) - (1 + y)*(1+t)*y_2prime_2) / (1 + y) ** 3
y_4prime_2 = (y_3prime_2*(-1*t - 2*t*y - t*y**2 - 1 - 2*y - y**2)
              + y_2prime_2*(6*t*y_1prime_2 + 6*t*y*y_1prime_2 - 3 - 6*y + 6*y_1prime_2 + 6*y*y_1prime_2 - 3*y**2)
              + (y_1prime_2**2)*(6 + 6*y) - (y_1prime_2**3)*(6 + 6*t)) / (1 + y) ** 4
a2 = 1
b2 = 2
n2 = 2
ic2 = 2
t2, y2 = taylor_two(y_1prime_2, y_2prime_2, var1, var2, a2, b2, n2, ic2)
print("Q2-b")
[print(i, end=' ') for i in t2]
print("")
[print(i, end=' ') for i in y2]
print("\n---------------------------------------------------")
t2, y2 = taylor_four(y_1prime_2, y_2prime_2, y_3prime_2, y_4prime_2, var1, var2, a2, b2, n2, ic2)
print("Q4-b")
[print(i, end=' ') for i in t2]
print("")
[print(i, end=' ') for i in y2]
print("\n---------------------------------------------------")
#################################################
y_1prime_3 = (2/t)*y + t**2 * E**t
y_2prime_3 = (2*(t*y_1prime_3 - y))/(t**2) + t**2 * E**t + 2*t * E**t
y_3prime_3 = (2*(t**2 * y_2prime_3 - 2*(t*y_1prime_3 - y)))/(t**3) + t**2 * E**t + 4*t * E**t + 2*E**t
y_4prime_3 = (2*(t**3 * y_3prime_3 - 3*(t**2 * y_2prime_3 - 2*(t*y_1prime_3 - y))))/(t**4)\
             + t**2 * E**t + 6*t * E**t + 6*E**t
a3 = 1
b3 = 2
n3 = 10
ic3 = 0
t3, y3 = taylor_two(y_1prime_3, y_2prime_3, var1, var2, a3, b3, n3, ic3)
print("Q9-a")
[print(i, end=' ') for i in t3]
print("")
[print(i, end=' ') for i in y3]
print("\n---------------------------------------------------")
t3, y3 = taylor_four(y_1prime_3, y_2prime_3, y_3prime_3, y_4prime_3, var1, var2, a3, b3, n3, ic3)
print("Q9-c")
[print(i, end=' ') for i in t3]
print("")
[print(i, end=' ') for i in y3]
print("\n---------------------------------------------------")
