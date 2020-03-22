import numpy as np
import math as mp
from matplotlib import pyplot as plt

tau = 10  # taux de compression
R = 0.038  # rayon du vilebrequin
D = 0.09  # dimaètre du piston/ cylindre
L = 0.117  # longueur de la bielle
Vc = mp.pi * (D ** 2) * R / 2  # volume dy cylindrée

Vmin = (1 / tau - 1) * Vc  # volume mort

Vmax = (tau / tau - 1) * Vc  # volume max

beta = L / R

gamma = 1.3

Qtot = (2800 * 10E3) * 1000 * Vmax #essence
theta_d = -mp.pi / 6
inter_theta = mp.pi / 4


def f(t, z):
    V = Vc * (((1 - mp.cos(t) + beta - mp.sqrt(beta ** 2 - (mp.sin(t)) ** 2)) / 2) + 1 / (tau - 1))

    dV = (Vc / 2) * (mp.sin(t) + ((mp.sin(t) * mp.cos(t)) / np.sqrt(beta ** 2 - (mp.sin(t)) ** 2)))

    dQ = (Qtot * mp.pi / 2 * inter_theta) * (mp.sin(mp.pi * (t - theta_d) / inter_theta))
    dp = ((-gamma * dV * z) / V) + ((gamma - 1) * dQ / V)
    return dp


def RG(Xstart, Xend, Ustart, n):
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

m_piston = - 0.01665 + 1.467E-6 * D**3 #masse piston diesel
m_bielle = 1.2 * m_piston  # source https://www.itterbeek.org/uploads/documents/BEMPr%C3%A9dimensionnement.pdf
RPM = 3000 #essence
w = 2*mp.pi /60 * RPM

def F_pied(X, Upressure):
    return ((mp.pi*D**2)/4)*Upressure - (m_piston*R * w**2 * np.cos(X))


def F_tete (X, Upressure):
    return  -((mp.pi*D**2)/4)*Upressure + (m_piston + m_bielle)*w**2 * np.cos(X)


plt.figure(figsize=(6, 9))

plt.subplot(3,1,1)
Xstart = -mp.pi; Xend =mp.pi
Ustart = 3
n = 1000
X, U = RG(Xstart, Xend, Ustart, n)
plt.title("Evolution de la pression")
plt.plot(X, U, '-r', linewidth=0.5)

plt.subplot(3,1,2)
F_pied = F_pied(X,U)
plt.title("Force sur le pied de la bielle")
plt.plot(X,F_pied)

plt.subplot(3,1,3)
F_tete = F_tete(X,U)
plt.title("Force sur la tete de la bielle")
plt.plot(X,F_tete)

plt.show()
