# 3 лаба пока не готова

import numpy
import matplotlib.pyplot as plotter


def f(x):
    return (numpy.pi * x - 3) / (x - 1) ** 2


def f2(x):
    return 2*(x-1) - 5*(numpy.pi*(x**2 + x - 2) - 9*x + 9) # 2⋅(x−1)−5⋅(π⋅(x2+x−2)−9⋅x+9)


def df(x):
    return 24*((x-1)**-7)*(numpy.pi*(x*(x+3) - 4)- 15*x + 15)


def Lagrange(arg, x, y):
    L = 0
    size = len(x)
    for i in range(size):
        l = 1
        for j in range(size):
            if i != j:
                l = l * (arg - x[j]) / (x[i] - x[j])
        L = L + l * y[i]
    return L


def Newton(x, X, Y, h):
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


a = 0
b = 0.3
h = (b - a) / 3

plotter.figure(1)
x = numpy.arange(a, b, 0.001)
plotter.title("График 4 производной")
plotter.xlabel("x")
plotter.ylabel("y")
plotter.grid()
plotter.plot(x, df(x))

X = numpy.array([a, a + h, a + 2 * h, a + 3 * h])
Y = f(X)
lagrange, newton, function = [], [], []

for i in range(len(x)):
    lagrange.append(Lagrange(x[i], X, Y))
    newton.append(Newton(x[i], X, Y, h))
    function.append(f(x[i]))

plotter.figure(2)
plotter.title("Полиномы")
plotter.xlabel("x")
plotter.grid()
plotter.plot(x, newton, label="Ньютон")
plotter.plot(x, lagrange, label="Лагранж")
plotter.plot(x, function, label="Функция")
plotter.legend()

plotter.figure(3)
plotter.title("Абслоютные погрешности")
plotter.xlabel("x")
plotter.ylabel("Rn(x)")
plotter.plot(x, abs(f(x) - newton), label="Ньютон")
plotter.plot(x, abs(f(x) - lagrange), label="Лагранж")
plotter.legend()
plotter.grid()
plotter.show()