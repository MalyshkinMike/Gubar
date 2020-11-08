# первая лаба
import matplotlib.pyplot as canvas
import numpy as np
from numpy import fabs


def func(x):
    return (x ** 3 - 0.9 * x ** 2 - 22.0 * x - 30.8)


def dichotomy(a, b, f, eps):
    i = 0
    while (fabs(b - a) > 2 * eps):
        i += 1
        x1 = (b + a) / 2
        y = f(x1)
        ymin = f(a)
        if ((y > 0 and ymin > 0) or (y < 0 and ymin < 0)):
            a = x1
        else:
            b = x1
    return x1, i



def gybrid(x_1, x0, f, eps):
    i = 0
    x1 = x0 - f(x0)(x0 - x_1)/(f(x0) - f(x_1))
    while True:
        i += 1
        if fabs(x1 - x0) <= eps:
            return x1, i

        x0 = x1
        x_1 = x0
        x1 = x0 - f(x0)(x0 - x_1)/(f(x0) - f(x_1))


gr_x_min = int(input("левая граница графика ►"))
gr_x_max = int(input("правая граница графика ►"))
gr_x = np.linspace(gr_x_min, gr_x_max, num=100)
gr_y = func(gr_x)
zero = [0] * 100
canvas.plot(gr_x, gr_y)
canvas.plot(gr_x, zero, color='black')  # ось х
canvas.xlabel('x')
canvas.ylabel('y')
canvas.grid()
canvas.show()

eps = 0.00001
left_border = float(input("левая граница начального приближения ► "))
right_border = float(input("правая граница начального приближения ► "))
dichotomy_res = dichotomy(left_border, right_border, func, eps)
print("Введите 2 точки (х) начального приближения для гибридного метода")
x_1 = float(input('\t1 ► '))
x0 = float(input('\t2 ► '))
gybrid_res = gybrid(x_1, x0, func, eps)
print("Метод дихотомии\n" + "\tРезультат ►" + str(dichotomy_res[0]) + "\n\tКоличество итераций ►" + str(
    dichotomy_res[1]) + "\n\tточность ►" + str(eps))
print("Гибридный метод\n" + "\tРезультат ►" + str(gybrid_res[0]) + "\n\tКоличество итераций ►" + str(
    gybrid_res[1]) + "\n\tточность ►" + str(eps))
