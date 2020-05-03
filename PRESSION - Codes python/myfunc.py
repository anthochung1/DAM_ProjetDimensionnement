import numpy as np
import math as mp
from matplotlib import pyplot as plt

C = 0.03  # 2 * R  # [m]
R = C / 2  # 0.038
tau = 10  # [-]
D = 0.05  # [m]
L = 0.35  # 0.117  # [m]
mpiston = 0.7  # [kg]
mbielle = 0.7  # [kg]
Q = 2800000  # [J/kg_inlet gas]

Vc = (np.pi * (D ** 2) * R) / 2  # volume dy cylindrée

beta = L / R

gamma = 1.3


def function(t, z, s, thetaC, deltaThetaC):
    if (t < thetaC % (2 * np.pi) or t > (thetaC + deltaThetaC) % (2 * np.pi)):
        Q_tot = 0
    else:
        Q_tot = (Q * ((s * (101325) * Vc) / (8.314 * 303.15))) * 0.029

    V = Vc * (((1 - np.cos(t) + beta - np.sqrt(beta ** 2 - (np.sin(t)) ** 2)) / 2) + 1 / (tau - 1))
    dV = (Vc / 2) * (np.sin(t) + (np.sin(t) * np.cos(t)) / np.sqrt(beta ** 2 - (np.sin(t)) ** 2))

    dQ = (Q_tot * np.pi / (2 * deltaThetaC)) * (np.sin(np.pi * (t - thetaC) / deltaThetaC))
    dp = (-gamma * dV * z + (gamma - 1) * dQ) / V

    return dp


def RG_function(Xstart, Xend, Ustart, n, s, thetaC, deltaThetaC):
    h = (Xend - Xstart) / n
    T = np.linspace(Xstart, Xend, n + 1)
    U = np.zeros(n + 1)
    U[0] = Ustart
    for i in range(n):
        K1 = function(T[i], U[i], s, thetaC, deltaThetaC)
        K2 = function(T[i] + h / 2, U[i] + K1 * h / 2, s, thetaC, deltaThetaC)
        K3 = function(T[i] + h / 2, U[i] + K2 * h / 2, s, thetaC, deltaThetaC)
        K4 = function(T[i] + h, U[i] + K3 * h, s, thetaC, deltaThetaC)
        U[i + 1] = U[i] + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
    return T, U


def F_tete_function(X, Upressure, rpm):
    w = (2 * np.pi / 60) * rpm
    return (-(np.pi * D ** 2) / 4) * Upressure + (mpiston + mbielle) * R * (w ** 2) * np.cos(X)


def F_pied_function(X, Upressure, rpm):
    w = (2 * np.pi / 60) * rpm
    return ((np.pi * D ** 2) / 4) * Upressure - (mpiston * R * (w ** 2) * np.cos(X))


def myfunc(rpm, s, thetaC, deltaThetaC):
    plt.figure(figsize=(6, 9))

    plt.subplot(3, 1, 1)
    Xstart = -np.pi
    Xend = 4 * np.pi
    Ustart = s * 101325
    n = 10000
    X, U = RG_function(Xstart, Xend, Ustart, n, s, -thetaC / 180 * np.pi, deltaThetaC / 180 * np.pi)

    plt.title("Evolution de la pression")
    plt.plot(X, U / (10 ** 5), '-r', linewidth=0.5)

    plt.subplot(3, 1, 2)
    F_pied = F_pied_function(X, U, rpm)
    plt.title("Force sur le pied de la bielle")
    plt.plot(X, F_pied)

    plt.subplot(3, 1, 3)
    F_tete = F_tete_function(X, U, rpm)

    plt.title("Force sur la tete de la bielle")
    plt.plot(X, F_tete)

    plt.show()

    F_crit = np.max(np.abs(F_tete) + np.abs(F_pied))
    sig = 450000000
    E = 200000000000

    a_xx = 1 / F_crit
    b_xx = -1 / (11 * sig)
    c_xx = -(0.5 * L) ** 2 / ((np.pi) ** 2 * E * (419 / 12))

    a_yy = 1 / F_crit
    b_yy = -1 / (11 * sig)
    c_yy = -(1 * L) ** 2 / ((np.pi) ** 2 * E * (131 / 12))

    coeff_x = [a_xx, 0, b_xx, 0, c_xx]
    coeff_y = [a_yy, 0, b_yy, 0, c_yy]

    root_array_x = np.roots(coeff_x)
    root_array_y = np.roots(coeff_y)
    root_x = [i.real for i in root_array_x if (i.imag == 0 and i >= 0)]
    root_y = [i.real for i in root_array_y if (i.imag == 0 and i >= 0)]

    return np.max(np.array([root_x, root_y]))


a = myfunc(1141, 1.9, 20, 65)
print(a)
