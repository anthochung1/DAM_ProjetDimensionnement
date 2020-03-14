import numpy as np
import math as mp

tau = 10  # taux de compression
R = 0.038  # rayon du vilebrequin
D = 0.09  # dimaètre du piston/ cylindre
L = 0.117  # longueur de la bielle
Vc = mp.pi * (D ** 2) * R / 2  # volume dy cylindrée

Vmin = (1 / tau - 1) * Vc  # volume mort

Vmax = (tau / tau - 1) * Vc  # volume max

beta = L / R

gamma = 1.3

Qtot = (2800 * 10E3) * 1000 * Vmax
theta_d = -mp.pi / 6
inter_theta = mp.pi / 4


def f(t, z):
    V = Vc * (((1 - mp.cos(t) + beta - mp.sqrt(beta ** 2 - (mp.sin(t)) ** 2)) / 2) + 1 / (tau - 1))

    dV = (Vc / 2) * (mp.sin(t) + ((mp.sin(t) * mp.cos(t)) / np.sqrt(beta ** 2 - (mp.sin(t)) ** 2)))

    dQ = (Qtot * mp.pi / 2 * inter_theta) * (mp.sin(mp.pi * (t - theta_d) / inter_theta))
    dp = ((-gamma * dV * z) / V) + ((gamma - 1) * dQ / V)
    return dp


def dim(Xstart, Xend, Ustart, n):
    h = (Xend - Xstart) / n
    T = np.linspace(Xstart, Xend, n + 1)
    U = np.zeros(n + 1)
    U[0] = Ustart
    for i in range(n):
        K1 = f(T[i], U[i])
        K2 = f(T[i] + h / 2, U[i] + K1 * h / 2)
        K3 = f(T[i] + h / 2, U[i] + K2 * h / 2)
        K4 = f(T[i] + h, U[i] + K3 * h)
        U[i + 1] = U[i] + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
    return T, U


"""def euler(Xstart, Xend, Ustart, n):
    T = np.linspace(Xstart, Xend, n + 1)
    U = np.zeros(n + 1)
    U[0] = Ustart
    for i in range(n):
        U[i + 1] = U[i] + f(T[i], U[i])
    return U"""

from matplotlib import pyplot as plt

plt.figure("Dimensionnement")
Xstart = -6.28; Xend =6
Ustart = 2*10**5
n = 1000

X, U = dim(Xstart, Xend, Ustart, n)
# X,U = dim( -4 * mp.pi,4 * mp.pi, 101325,10000)

"""T = np.linspace(Xstart, Xend, n + 1)
Y = euler(Xstart, Xend, Ustart, n)"""

plt.plot(X, U, '-r', linewidth=0.5)
"""plt.plot(T, Y, '-g', linewidth=0.5)"""
plt.show()
print(U)
