# 5 готова, но не маскирована под оригинальность

import numpy
import matplotlib.pyplot as canvas
from numpy import log as ln


def func(x):
    return (x ** 2) * ln(x)


def sec_diff(x):  # 2xln(x) + x**2/x = 2xln(x) + x
    # 2ln(x) + 2 + 1
    return 2 * ln(x) + 3


def trap_method(x, h, a, b):
    size = len(x)
    res = (func(a) + func(b)) / 2
    for i in range(1, size):
        res += func(x[i])
    return res * h


def simpson_method(x, h, a, b):
    size = int(len(x) / 2)
    res = (func(a) + 4 * func(a + h) + func(b))
    for i in range(1, size):
        res += 2 * func(x[2 * i]) + 4 * func(x[2 * i + 1])
    return res * h / 3


a = 1
b = 2
d = numpy.arange(a, b, 0.001)
canvas.figure(1)
canvas.title("2я производная")
canvas.xlabel("Х")
canvas.ylabel("Y")
canvas.grid()
canvas.plot(d, sec_diff(d))
x_maximum = 2.0  # максимум(для значения 2й производной ) - по графику
m = sec_diff(x_maximum)  # Максимальное значение 2й производной на промежутке
# eps = 0.00001
# h = numpy.sqrt(eps*12/((b-a)*m)) # h = 0.005230482293837083
# n = (b-a)/h #n = 192 (на самом деле 191.186958...., но 192 кратно 4
h = (b - a) / 192
x1 = numpy.arange(a, b, h)
x2 = numpy.arange(a, b, 2 * h)
trap_res1 = trap_method(x1, h, a, b)
trap_res2 = trap_method(x2, 2 * h, a, b)
xs1 = numpy.arange(a, b, h)
xs2 = numpy.arange(a, b, 2 * h)
sim_res1 = simpson_method(xs1, h, a, b)
sim_res2 = simpson_method(xs2, 2 * h, a, b)
exact_res = 1.070614703715409714001507879444  # точный результат
print("Для метода трапеций:")
print("Результат расчета с шагом h: " + str(trap_res1))
print("Результат расчета с шагом 2h: " + str(trap_res2))
print("Сравнение с точным результатом: " + str(abs(exact_res - trap_res1)))
print("Погрешность по правилу Рунге: " + str(abs(trap_res2 - trap_res1) / 3))
print("Для метода Симпсона:")
print("Результат расчета с шагом h: " + str(sim_res1))
print("Результат расчета с шагом 2h: " + str(sim_res2))
print("Сравнение с точным результатом: " + str(abs(exact_res - sim_res1)))
print("Погрешность по правилу Рунге: " + str(abs(sim_res2 - sim_res1) / 15))
canvas.show()
