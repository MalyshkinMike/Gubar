import numpy as np
import matplotlib.pyplot as canvas


def f(x):
    return (np.pi * x - 3) / (x - 1) ** 2


def f2(x):
    return 2 * (x - 1) - 5 * (np.pi * (x ** 2 + x - 2) - 9 * x + 9)  # 2⋅(x−1)−5⋅(π⋅(x2+x−2)−9⋅x+9)


def f4(x):
    return 24 * ((x - 1) ** -7) * (np.pi * (x * (x + 3) - 4) - 15 * x + 15)


def linear_S(dots, x, y):
    h = x[1] - x[0]
    S = np.array([y[0],
                  y[1],
                  y[2],
                  (y[1] - y[0]) / h,
                  (y[2] - y[1]) / h,
                  (y[3] - y[2]) / h
                  ])
    res = []
    for i in range(len(dots)):
        if x[0] <= dots[i] <= x[1]:
            res.append(S[0] + S[3] * (dots[i] - x[0]))
        elif x[1] <= dots[i] <= x[2]:
            res.append(S[1] + S[4] * (dots[i] - x[1]))
        else:
            res.append(S[2] + S[5] * (dots[i] - x[2]))
    return res


def parabol_S(dots, x, y):
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


def cube_S(dots, x, y):
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


def lagrange_P(arg, x, y):
    L = 0
    size = len(x)
    for i in range(size):
        l = 1
        for j in range(size):
            if i != j:
                l = l * (arg - x[j]) / (x[i] - x[j])
        L = L + l * y[i]
    return L


def newton_P(x, X, Y, h):
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
b = 0.43
x = np.arange(a, b, 0.001)
# h = 0.16
h = (b - a) / 3
X = np.array([a, a + h, a + 2 * h, a + 3 * h])
Y = f(X)
lagrange, newton, func = [], [], []
for i in range(len(x)):
    lagrange.append(lagrange_P(x[i], X, Y))
    newton.append(newton_P(x[i], X, Y, h))
    func.append(f(x[i]))
linear_spline = linear_S(x, X, Y)
parabol_spline = parabol_S(x, X, Y)
cube_spline = cube_S(x, X, Y)
canvas.figure(4)
canvas.title("Сплайны")
canvas.xlabel("x")
canvas.ylabel("y")
canvas.grid()
canvas.plot(x, func, label="Функция")
canvas.plot(x, linear_spline, label="L Сплайн")
canvas.plot(x, parabol_spline, label="P Сплайн")
canvas.plot(x, cube_spline, label="C Сплайн")
canvas.legend()
canvas.figure(5)
canvas.title("Абслоютные погрешности")
canvas.xlabel("x")
canvas.ylabel("Rn(x)")
canvas.plot(x, abs(f(x) - newton), label="Ньютон")
canvas.plot(x, abs(f(x) - lagrange), label="Лагранж")
canvas.plot(x, abs(f(x) - linear_spline), label="L Сплайн")
canvas.plot(x, abs(f(x) - parabol_spline), label="P Сплайн")
canvas.plot(x, abs(f(x) - cube_spline), label="C Сплайн")
canvas.legend()
canvas.grid()
canvas.show()
