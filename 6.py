# 6 лаба тоже рабочая, но не "заоригиналеная".
import numpy as np
import matplotlib.pyplot as plt
a = 0
b = 1
n = 8
x0 = 0
y0 = 1
EPS = 0.001

def f(x, y):
    return 2*x**3*y**3 - 2*x*y
def correct(x):
    return np.sqrt(2/(2*x**2 + np.exp(2*x**2) + 1))
def delta_f(h, x, y):
    ph1 = f(x, y)
    ph2 = f(x + h / 2, y + h * ph1 / 2)
    ph3 = f(x + h / 2, y + h * ph2 / 2)
    ph4 = f(x + h, y + h * ph3)
    return h * 1 / 6 * (ph1 + 2 * ph2 + 2 * ph3 + ph4)
def find_h(y0, x0, h):
    y1 = y0 + delta_f(h, x0, y0)
    y2 = y1 + delta_f(h, x0 + h, y1)
    y2kr = y0 + delta_f(2 * h, x0, y0)
    return abs(y2 - y2kr)
def runge(h):
    xi = np.linspace(a, b, num=n + 1)
    yi = [y0]
    for i in range(len(xi)):
        yi.append(yi[i] + delta_f(h, xi[i], yi[i]))
    yi.pop(-1)
    return yi
def eiler(h):
    yi = [y0]
    xi = np.linspace(a, b, num=n + 1)
    for i in range(len(xi)):
        yi.append(yi[i] + h * f(xi[i], yi[i]))
    yi.pop(-1)
    return yi
def double_step_runge(h, runge_y):
    xi = np.linspace(a, b, num=n + 1)
    yi = [y0]

    dyi = []
    j = 0
    for i in range(0, len(xi), 2):
        yi.append(delta_f(2*h, xi[i], yi[j]))
        dyi.append(runge_y[i] - yi[j])
        j += 1
    yi.pop(-1)
    return yi, dyi
def doube_step_eiler(h, eil_y):
    xi = np.linspace(a, b, num=n + 1)
    yi = [y0]
    dyi = []
    j = 0
    for i in range(0, len(xi), 2):
        yi.append(yi[j] + 2*h * f(xi[i], yi[j]))
        dyi.append(eil_y[i] - yi[j])
        j += 1
    yi.pop(-1)
    return yi, dyi
def compare(rung_y, eil_y):
    x = np.linspace(a, b, num=n + 1)
    runge_y_wave, runge_delta_y = double_step_runge(h, rung_y)
    eiler_y_wave, eiler_delta_y = doube_step_eiler(h, eil_y)
    print("Рунге-Кутта")
    print(f" x\t\t yi\t\tyi2\t\tdy")
    for i in range(n):
        if ( i%2 == 0):
            print("%8.5f %8.5f %8.5f %8.5f" % (x[i], rung_y[i],
            runge_y_wave[i // 2], runge_delta_y[i // 2]))
        else:
            print("%8.5f %8.5f" % (x[i], rung_y[i]))
    print("Эйлер")
    print(f" x\t\t yi\t\tyi2\t\tdy")
    for i in range(n):
        if (i % 2 == 0):
            print("%8.5f %8.5f %8.5f %8.5f" % (x[i], eiler_y[i],
            eiler_y_wave[i // 2], eiler_delta_y[i // 2]))
        else:
            print("%8.5f %8.5f" % (x[i], eiler_y[i]))
def R(function_y):
    max = 0
    x = 0
    for i in range(n):
        if max < abs(function_y[i] - correct(x + i * 0.2)):
            max = abs(function_y[i] - correct(x + i * 0.2))
    return max
h = 0.1
while True:
    e = find_h(1.1, 0, h)
    while e < 10 ** -4:
        print("e = " + str(e))
        print(f"h = {h}")
        print("n = " + str(0.8 / h))
        h *= 2
        e = find_h(1.1, 0, h)
    print("e = " + str(e))
    print(f"h = {h}")
    print("n = " + str(0.8 / h))
    n = int((b - a) / h)
    if (n%2 == 0):
        if (e > EPS):
            print(f"последний полученный h является неточным. Возьмем предыдущий, (h = {h/2} )")
            h = h / 2
        else:
            print(f"h = {h}")
        break
    else:
        n += 1
        h = (b - a) / n
runge_y = runge(h)
eiler_y = eiler(h)
print("Максимальная погрешность для Рунге:" + str(R(runge_y)))
print("Максимальная погрешность для Эйлера:" + str(R(eiler_y)))
print()
print("Сравнение решений в узловых точках: ")
compare(runge_y, eiler_y)
x = np.linspace(a, b, num=n + 1)
plt.figure(1)
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.plot(x, runge_y, label="Рунге-Кутта")
plt.plot(x, eiler_y, label="Эйлера")
plt.plot(x, correct(x), label="Точное")
plt.legend()
plt.show()