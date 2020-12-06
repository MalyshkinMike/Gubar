# 6 лаба
import numpy as np
import matplotlib.pyplot as canvas

a = 0
b = 1
n = 8
x0 = 0
y0 = 1
eps = 0.001


def f(x, y):
    return 2 * x ** 3 * y ** 3 - 2 * x * y


def exact(x):
    return np.sqrt(2 / (2 * x ** 2 + np.exp(2 * x ** 2) + 1))


def Runge_Kutta_helper(h, x, y):
    ph1 = f(x, y)
    ph2 = f(x + h / 2, y + h * ph1 / 2)
    ph3 = f(x + h / 2, y + h * ph2 / 2)
    ph4 = f(x + h, y + h * ph3)
    return h * 1 / 6 * (ph1 + 2 * ph2 + 2 * ph3 + ph4)


def h_helper(y0, x0, h):
    y1 = y0 + Runge_Kutta_helper(h, x0, y0)
    y2 = y1 + Runge_Kutta_helper(h, x0 + h, y1)
    y2kr = y0 + Runge_Kutta_helper(2 * h, x0, y0)
    return abs(y2 - y2kr)


def Runge_Kutta(h):
    xi = np.linspace(a, b, num=n + 1)
    yi = [y0]
    for i in range(len(xi)):
        yi.append(yi[i] + Runge_Kutta_helper(h, xi[i], yi[i]))
    yi.pop(-1)
    return yi


def Eiler(h):
    yi = [y0]
    xi = np.linspace(a, b, num=n + 1)
    for i in range(len(xi)):
        yi.append(yi[i] + h * f(xi[i], yi[i]))
    yi.pop(-1)
    return yi


def double_step_Runge(h, runge):
    xi = np.linspace(a, b, num=n + 1)
    yi = [y0]

    y_difference = []
    j = 0
    for i in range(0, len(xi), 2):
        yi.append(Runge_Kutta_helper(2 * h, xi[i], yi[j]))
        y_difference.append(runge[i] - yi[j])
        j += 1
    yi.pop(-1)
    return yi, y_difference


def doube_step_Eiler(h, eiler):
    xi = np.linspace(a, b, num=n + 1)
    yi = [y0]
    y_difference = []
    j = 0
    for i in range(0, len(xi), 2):
        yi.append(yi[j] + 2 * h * f(xi[i], yi[j]))
        y_difference.append(eiler[i] - yi[j])
        j += 1
    yi.pop(-1)
    return yi, y_difference


def comparing(runge, eiler):
    x = np.linspace(a, b, num=n + 1)
    runge_y_wave, runge_delta_y = double_step_Runge(h, runge)
    eiler_y_wave, eiler_delta_y = doube_step_Eiler(h, eiler)
    print("Рунге-Кутта")
    print(f" x\t\t yi\t\tyi2\t\tdy")
    for i in range(n):
        if (i % 2 == 0):
            print("%8.5f %8.5f %8.5f %8.5f" % (x[i], runge[i],
                                               runge_y_wave[i // 2], runge_delta_y[i // 2]))
        else:
            print("%8.5f %8.5f" % (x[i], runge[i]))
    print("Эйлер")
    print(f" x\t\t yi\t\tyi2\t\tdy")
    for i in range(n):
        if (i % 2 == 0):
            print("%8.5f %8.5f %8.5f %8.5f" % (x[i], eiler[i],
                                               eiler_y_wave[i // 2], eiler_delta_y[i // 2]))
        else:
            print("%8.5f %8.5f" % (x[i], eiler[i]))


def R(function_y):
    max = 0
    x = 0
    for i in range(n):
        if max < abs(function_y[i] - exact(x + i * 0.2)):
            max = abs(function_y[i] - exact(x + i * 0.2))
    return max


h = 0.1
while True:
    e = h_helper(1.1, 0, h)
    while e < 10 ** -4:
        print("e = " + str(e))
        print(f"h = {h}")
        print("n = " + str(0.8 / h))
        h *= 2
        e = h_helper(1.1, 0, h)
    print("e = " + str(e))
    print(f"h = {h}")
    print("n = " + str(0.8 / h))
    n = int((b - a) / h)
    if (n % 2 == 0):
        if (e > eps):
            print(f"последний полученный h является неточным. Возьмем предыдущий, (h = {h / 2} )")
            h = h / 2
        else:
            print(f"h = {h}")
        break
    else:
        n += 1
        h = (b - a) / n
runge_res = Runge_Kutta(h)
eiler_res = Eiler(h)
print("Максимальная погрешность для Рунге:" + str(R(runge_res)))
print("Максимальная погрешность для Эйлера:" + str(R(eiler_res)))
print()
print("Сравнение решений в узловых точках: ")
comparing(runge_res, eiler_res)
x = np.linspace(a, b, num=n + 1)

canvas.xlabel("x")
canvas.ylabel("y")
canvas.grid()
canvas.plot(x, eiler_res, label="Эйлера")
canvas.plot(x, runge_res, label="Runge")
canvas.plot(x, exact(x), label="Точное")
canvas.legend()
canvas.show()