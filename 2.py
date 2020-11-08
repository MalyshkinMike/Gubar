from numpy import cos, sin, fabs
from matplotlib import pyplot as canvas
import numpy as np

# 2 лаба
eq1 = lambda x, y: cos(x + 0.5) - 2 - y  # eq - equation - уравнение
eq2 = lambda x, y: sin(y) - 2 * x - 1
equations = lambda x, y: np.array([[eq1(x, y)], [eq2(x, y)]])


def graph():
    f1 = lambda x: cos(x + 0.5) - 2
    f2 = lambda y: (sin(y) - 1) / 2
    x = np.linspace(-2, 0, 200)
    canvas.xlabel('X')
    canvas.ylabel('Y')
    canvas.grid()

    canvas.plot(x, f1(x), label='cos(x+0.5) - y = 2')
    canvas.plot(f2(x), x, label='sin(y) - 2x = 1')
    canvas.legend()


J = np.array([
    [lambda x, y: -sin(x + 0.5), lambda x, y: -1],
    [lambda x, y: -2, lambda x, y: cos(y)]
])


def J_inv(X, Y):
    determinant = lambda x, y: -2 - sin(x + 0.5) * cos(y)
    return np.array([
        [(lambda x, y: cos(y) / determinant(x, y))(X, Y), (lambda x, y: 1 / determinant(x, y))(X, Y)],
        [(lambda x, y: 2 / determinant(x, y))(X, Y), (lambda x, y: -sin(x + 0.5) / determinant(x, y))(X, Y)]
    ])


def Newton_method(x0, y0, eps):
    i = 0
    [[x1], [y1]] = np.array([[x0], [y0]]) - J_inv(x0, y0).dot(equations(x0, y0))
    pre_x, pre_y = x0, y0
    while True:
        if (fabs(x1 - x0) <= eps and fabs(y1 - y0) <= eps):
            return x1, y1, i
        if (fabs(x1 - x0) >= fabs(pre_x - x0) or fabs(y1 - y0) >= fabs(pre_y - y0)) and i != 0:
            return np.nan, np.nan, np.nan
        pre_x = x0
        pre_y = y0
        x0 = x1
        y0 = y1
        [[x1], [y1]] = np.array([[x0], [y0]]) - J_inv(x0, y0).dot(equations(x0, y0))
        i += 1


def modified_Newton_method(x0, y0, eps):
    i = 0
    [[x1], [y1]] = np.array([[x0], [y0]]) - J_inv(x0, y0).dot(equations(x0, y0))
    x00 = x0
    y00 = y0
    pre_x, pre_y = x0, y0
    while True:
        if (fabs(x1 - x00) <= eps and fabs(y1 - y00) <= eps):
            return x1, y1, i
        if (fabs(x1 - x00) >= fabs(pre_x - x00) or fabs(y1 - y00) >= fabs(pre_y - y00)) and i != 0:
            return np.nan, np.nan, np.nan
        pre_x, pre_y = x00, y00
        x00 = x1
        y00 = y1
        [[x1], [y1]] = np.array([[x00], [y00]]) - J_inv(x0, y0).dot(equations(x00, y00))
        i += 1


graph()
canvas.show()
eps = 0.00001
print('cos(x+0.5) - y = 2')
print('sin(y) - 2x = 1')
print('Начальное приближение ► ')
x = float(input('\tВведите х ► '))
y = float(input('\tВведите y ► '))
newt_x, newt_y, newt_i = Newton_method(x, y, eps)
mod_newt_x, mod_newt_y, mod_newt_i = modified_Newton_method(x, y, eps)
print('Метод Ньютона')
print('\tx = ' + str(newt_x) + '\n\ty = ' + str(newt_y) + '\n\tколичество итераций (i) = ' + str(newt_i))
print('Метод Ньютона (модифицированный)')
print('\tx = ' + str(mod_newt_x) + '\n\ty = ' + str(mod_newt_y) + '\n\tколичество итераций (i) = ' + str(mod_newt_i))
