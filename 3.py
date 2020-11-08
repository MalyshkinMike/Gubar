# 3 лаба

import numpy as np
from numpy import fabs
import matplotlib.pyplot as canvas


def f(x):
    return (np.pi * x - 3) / (x - 1) ** 2


def f2(x):
    return 2 * (x - 1) - 5 * (np.pi * (x ** 2 + x - 2) - 9 * x + 9)  # 2⋅(x−1)−5⋅(π⋅(x2+x−2)−9⋅x+9)


def f4(x):
    return 24 * ((x - 1) ** -7) * (np.pi * (x * (x + 3) - 4) - 15 * x + 15)


def f4abs(x):
    return fabs(f4(x))


def lagrange_P(arg, x, y):
    Lagrange = 0
    size = len(x)
    for i in range(size):
        l = 1
        for j in range(size):
            if i != j:
                l = l * (arg - x[j]) / (x[i] - x[j])
        Lagrange = Lagrange + l * y[i]
    return Lagrange


def newton_P(x, X, Y, h):
    Newton = Y[0]
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
        Newton = Newton + temp
    return Newton


def func_max(func, x1, x2, eps):
    while fabs(x1 - x2) >= eps:
        y1 = fabs(func(((x1 + x2) / 2) - eps))
        y2 = fabs(func(((x1 + x2) / 2) + eps))
        if y1 > y2:
            x2 = x1 + ((x2 - x1) / 2) - eps
        else:
            x1 = x1 + ((x2 - x1) / 2) + eps
    else:
        return (x1 + x2) / 2


a = 0
eps = 0.001
b = 0.9

'''def find_err(a_arg, b_arg):
    xx = np.linspace(a_arg, b_arg, num=4)

    def w_4(argument):
        return (argument - xx[0]) * (argument - xx[1]) * (argument - xx[2]) * (argument - xx[3])

    return find_max(w_4, a_arg, b_arg, eps) * find_max(f4abs, a_arg, b_arg, eps) / 24


e = find_err(a, b)

while e >= eps:
    b -= 0.01
    e = find_err(a, b)
print(b)
путем сужения отрезка получено b = 0.43
'''

b = 0.43
h = (b - a) / 3
canvas.figure(1)
x = np.arange(a, b, 0.001)
canvas.title("График 4 производной")
canvas.xlabel("x")
canvas.ylabel("y")
canvas.grid()
canvas.plot(x, f4(x))

X = np.array([a, a + h, a + 2 * h, a + 3 * h])
Y = f(X)
lagrange, newton, function = [], [], []

for i in range(len(x)):
    lagrange.append(lagrange_P(x[i], X, Y))
    newton.append(newton_P(x[i], X, Y, h))
    function.append(f(x[i]))

canvas.figure(2)
canvas.title("Полиномы")
canvas.xlabel("x")
canvas.grid()
canvas.plot(x, newton, label="Ньютон")
canvas.plot(x, lagrange, label="Лагранж")
canvas.plot(x, function, label="Функция")
canvas.legend()

canvas.figure(3)
canvas.title("Абслоютные погрешности")
canvas.xlabel("x")
canvas.ylabel("Rn(x)")
canvas.plot(x, abs(f(x) - newton), label="Ньютон")
canvas.plot(x, abs(f(x) - lagrange), label="Лагранж")
canvas.legend()
canvas.grid()
canvas.show()
