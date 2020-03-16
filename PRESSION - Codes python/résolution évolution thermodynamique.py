import numpy as np
import math as m
import matplotlib.pyplot as plt


def euler_explicite(f, init, end, n, C):
    """
    :param f: fonction correspondant à u' = f(x,u(x))
    :param init: tuple contenant (float abscisse_de_départ, float ordonée_de_départ)
    :param end: float correspondant à l'abscisse finale pour laquelle sera évaluée u.
    :param n: nombre entier de pas que va effectuer la méthode.
    :return: un tuple de deux arrays de n éléments (array_des_abscisses, array_des_ordonnées)
    """

    X = np.linspace(init[0], end, n)
    print(X)
    dx = X[1] - X[0]
    Y = np.zeros(n)
    Y[0] = init[1]
    for i in range(len(X) - 1):
        Y[i + 1] = Y[i] + f(X[i], Y[i], C) * dx
    return X, Y


def evol(t, p, C):
    """
    :param t: angle du vilbrequin (en radian)
    :param p: pression correspondante à l'angle t du vilbrequin en Pascal
    :param C: Dictionnaire des paramètres de la fonction
    :return: dp/dt l'évolution de la pression en fonction de l'angle t du vilbrequin en Pascal
    """
    V = C["Vc"] * (((1 - np.cos(t) + C["beta"] - np.sqrt(C["beta"] ** 2 - (np.sin(t)) ** 2)) / 2) + 1 / (C["tau"] - 1))
    dV = C["Vc"] * (np.sin(t) + (np.sin(t) * np.cos(t)) / np.sqrt(C["beta"] ** 2 - (np.sin(t)) ** 2))

    dQ = (C["Qtot"] * np.pi / 2 * C["Dtheta"]) * (np.sin(np.pi * (t - C["theta_d"]) / C["Dtheta"]))

    return -C["gamma"] * (p / V) * dV + ((C["gamma"] - 1) / V) * dQ


# Myparameters correspond aux différents paramètres du système.

Myparameters = {"Vc": 0.01, "beta": 3.5, "tau": 10, "Dtheta": np.pi / 3, "theta_d": np.pi / 6, "gamma": 1.3,
                "Qtot": 2800}
x, y = euler_explicite(evol, (-m.pi, 101325), m.pi, 10000, Myparameters)
plt.plot(x, y / 101325, "r-")
plt.show()
