# 5 лаба
import numpy as np
import matplotlib.pyplot as canvas
from numpy import log as ln


def f(x):
    return (x ** 2) * ln(x)


def f2(x):
    return 2 * ln(x) + 3


def trapecy(x, h, a, b):
    size = len(x)
    res = (f(a) + f(b)) / 2
    for i in range(1, size):
        res += f(x[i])
    return res * h


def simpson(x, h, a, b):
    size = int(len(x) / 2)
    res = (f(a) + 4 * f(a + h) + f(b))
    for i in range(1, size):
        res += 2 * f(x[2 * i]) + 4 * f(x[2 * i + 1])
    return res * h / 3


a = 1
b = 2
d = np.arange(a, b, 0.001)
canvas.figure(1)
canvas.title("2я производная")
canvas.xlabel("Х")
canvas.ylabel("Y")
canvas.grid()
canvas.plot(d, f2(d))
x_m = 2.0  # максимум 2й производной. Получен графическим методом
m = f2(x_m)  # максимальное значение 2й производной
# eps = 0.00001
# h = numpy.sqrt(eps*12/((b-a)*m)) # h = 0.005230482293837083
# n = (b-a)/h #n = 192 (на самом деле 191.186958...., но берем 192 т.к. кратно 4) 192/4 = 48
h = (b - a) / 192
xh = np.arange(a, b, h)
x2h = np.arange(a, b, 2 * h)
trapecy1 = trapecy(xh, h, a, b)
trapecy2 = trapecy(x2h, 2 * h, a, b)
simpson1 = simpson(xh, h, a, b)
simpson2 = simpson(x2h, 2 * h, a, b)
exact = 1.070614703715409714001507879444  # точный результат
print("Метод трапеций:")
print("\tРезультат с шагом h ► " + str(trapecy1) + "\n\t\tс шагом 2h ► " + str(trapecy2))
print("\tСравнение с точным ► " + str(abs(exact - trapecy1)))
print("\tПогрешность по Рунге ► " + str(abs(trapecy2 - trapecy1) / 3))
print("Метода Симпсона:")
print("\tРезультат с шагом h ► " + str(simpson1) + "\n\t\tс шагом 2h ► " + str(simpson2))
print("\tСравнение с точным ► " + str(abs(exact - simpson1)))
print("\tПогрешность по Рунге ► " + str(abs(simpson2 - simpson1) / 15))
canvas.show()
