# 4 лаба пока мой вариант, после 3 обновлю.

import numpy as np
import matplotlib.pyplot as plt
def func(arg):
    return arg ** 2 * np.pi * np.exp(5 * arg) + 1 / arg
def f2(arg):
    return 2 * np.pi * np.exp(5 * arg) + 20 * arg * np.pi * np.exp(5 * arg) + 25 * np.pi * arg ** 2 * np.exp(5 * arg) + 2 / arg ** 3
def dif_4(arg):
    return 625 * np.pi * arg ** 2 * np.exp(5 * arg) + 1000 * np.pi * arg * np.exp(5 * arg) + 300 * np.pi * np.exp(5 * arg) + 24 / (arg ** 5)
def lin_spline(dots, x, y):
    h = x[1] - x[0]
    S = np.array([y[0],
    y[1],
    y[2],
    (y[1] - y[0]) / h,
    (y[2] - y[1]) / h,
    (y[3] - y[2]) / h])
    res = []
    for i in range(len(dots)):
        if x[0] <= dots[i] <= x[1]:
            res.append(S[0] + S[3] * (dots[i] - x[0]))
        elif x[1] <= dots[i] <= x[2]:
            res.append(S[1] + S[4] * (dots[i] - x[1]))
        else:
            res.append(S[2] + S[5] * (dots[i] - x[2]))
    return res

def par_spline(dots, x, y):
    a = 2
# a1,a2,a3, b1,b2, b3,c1, c2, c3
    A = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, h, 0, 0, h * h, 0, 0],
    [0, 1, 0, 0, h, 0, 0, h * h, 0],
    [0, 0, 1, 0, 0, h, 0, 0, h * h],
    [0, 0, 0, 1, -1, 0, 2 * h, 0, 0],
    [0, 0, 0, 0, 1, -1, 0, 2 * h, 0],
    [0, 0, 0, 0, 0, 0, 2, 0, 0]])
    B = np.array([[y[0]],
    [y[1]],
    [y[2]],
    [y[1]],
    [y[2]],
    [y[3]],
    [0],
    [0],
    [f2(a)]])
    S = np.linalg.solve(A, B)
    S = S.ravel()
    res = []

    for i in range(len(dots)):
        if x[0] <= dots[i] <= x[1]:
            t = dots[i] - x[0]
            res.append(S[0] + S[3] * t + S[6] * t * t)
        elif x[1] <= dots[i] <= x[2]:
            t = dots[i] - x[1]
            res.append(S[1] + S[4] * t + S[7] * t * t)
        else:
            t = dots[i] - x[2]
            res.append(S[2] + S[5] * t + S[8] * t * t)
    return res

def cub_spline(dots, x, y):
    # a1,a2,a3,b1,b2, b3,c1, c2, c3, d1, d2, d3
    A = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, h, 0, 0, h * h, 0, 0, h ** 3, 0, 0],
    [0, 1, 0, 0, h, 0, 0, h * h, 0, 0, h ** 3, 0],
    [0, 0, 1, 0, 0, h, 0, 0, h * h, 0, 0, h ** 3],
    [0, 0, 0, 1, -1, 0, 2 * h, 0, 0, 3 * h * h, 0, 0],
    [0, 0, 0, 0, 1, -1, 0, 2 * h, 0, 0, 3 * h * h, 0],
    [0, 0, 0, 0, 0, 0, 2, -2, 0, 6 * h, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 6 * h, 0],
    [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 6 * h]])
    B = np.array([[y[0]],
    [y[1]],
    [y[2]],
    [y[1]],
    [y[2]],
    [y[3]],
    [0],
    [0],
    [0],
    [0],
    [f2(a)],
    [f2(b)]])
    S = np.linalg.solve(A, B)
    S = S.ravel()
    res = []
    for i in range(len(dots)):
        if x[0] <= dots[i] <= x[1]:
            t = dots[i] - x[0]
            res.append(S[0] + S[3] * t + S[6] * t ** 2 + S[9] * t ** 3)
        elif x[1] <= dots[i] <= x[2]:
            t = dots[i] - x[1]
            res.append(S[1] + S[4] * t + S[7] * t ** 2 + S[10] * t ** 3)
        else:
            t = dots[i] - x[2]
            res.append(S[2] + S[5] * t + S[8] * t ** 2 + S[11] * t ** 3)
    return res
def lagrange(arg, x, y):
    L = 0
    size = len(x)
    for i in range(size):
        l = 1
        for j in range(size):
            if i != j:
                l = l * (arg - x[j]) / (x[i] - x[j])
                L = L + l * y[i]

    return L

def newton(x, X, Y, h):
    res = Y[0]
    t = (x - X[0]) / h
    diff = [0] * 4
    for i in range(len(Y)):
        diff[i] = Y[i]
        for i in range(len(X) - 1):
            for j in range(len(X) - i - 1):
                diff[j] = diff[j + 1] - diff[j]
            temp = diff[0]
            for j in range(0, i + 1):
                temp *= (t - j)
                temp /= (j + 1)
            res = res + temp
    return res
a = 2
b = 2.5
x = np.arange(a, b, 0.001)
# h = 0.16
h = (b - a) / 3
X = np.array([a, a + h, a + 2 * h, a + 3 * h])
Y = func(X)
lagr, newt, function = [], [], []
for i in range(len(x)):
    lagr.append(lagrange(x[i], X, Y))
    newt.append(newton(x[i], X, Y, h))
    function.append(func(x[i]))
l_spline = lin_spline(x, X, Y)
p_spline = par_spline(x, X, Y)
c_spline = cub_spline(x, X, Y)
plt.figure(4)
plt.title("Сплайны")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.plot(x, function, label="Функция")
plt.plot(x, l_spline, label="L Сплайн")
plt.plot(x, p_spline, label="P Сплайн")
plt.plot(x, c_spline, label="C Сплайн")
plt.legend()
plt.figure(5)
plt.title("Абслоютные погрешности")
plt.xlabel("x")
plt.ylabel("Rn(x)")
plt.plot(x, abs(func(x) - newt), label="Ньютон")
plt.plot(x, abs(func(x) - lagr), label="Лагранж")
plt.plot(x, abs(func(x) - l_spline), label="L Сплайн")
plt.plot(x, abs(func(x) - p_spline), label="P Сплайн")
plt.plot(x, abs(func(x) - c_spline), label="C Сплайн")
plt.legend()
plt.grid()
plt.show()