import matplotlib.pyplot as plt
import numpy as np

### OBJECTIVE ###
"""We want to find out the evolution of the pressure in an engine cylinder 
on the basis of the following ordinary differential equation :
        
        dP/dT = -g (P/V) (dV/dT) + (g-1) 1/V dQ/dT
        
        where:
         * P : pressure in the cylinder
         * T : angle between the crank and the cylinder axis
         * V : volume occupied by the fluid
         * Q : Gain of heat
         
"""

### SOLVING METHOD ###
"""To solve the above presented ODE, we here use the Runge-Kutta-Fehlberg resolution method."""


## Make RKF45 as general as possible (whatever the function we choose) ##
def feval(funcName, *args):
    return eval(funcName)(*args)


## Implementation of the Runge-Kutta-Fehlberg method ##
def RKF45(func, y_init, x_range, h):
    # m = len(yinit)

    n = int((x_range[-1] - x_range[0]) / h)

    # Containers for solutions
    x = np.zeros(n)
    x[0] = x_range[0]

    y = np.zeros(n)
    y[0] = y_init

    # Find value at each discrete value
    for i in range(n - 1):
        k1 = feval(func, x[i], y[i])

        k2 = feval(func, x[i] + h / 5, y[i] + k1 * (h / 5))

        k3 = feval(func, x[i] + (3 * h / 10), y[i] + k1 * (3 * h / 40) + k2 * (9 * h / 40))

        k4 = feval(func, x[i] + (3 * h / 5), y[i] + k1 * (3 * h / 10) - k2 * (9 * h / 10) + k3 * (6 * h / 5))

        k5 = feval(func, x[i] + h,
                   y[i] - k1 * (11 * h / 54) + k2 * (5 * h / 2) - k3 * (70 * h / 27) + k4 * (35 * h / 27))

        k6 = feval(func, x[i] + (7 * h / 8), y[i] + k1 * (1631 * h / 55296) + k2 * (175 * h / 512)
                   + k3 * (575 * h / 13824) + k4 * (44275 * h / 110592) + k5 * (253 * h / 4096))

        y[i + 1] = y[i] + h * (37 * k1 / 378 + 250 * k3 / 621 + 125 * k4 / 594 + 512 * k6 / 1771)
        x[i + 1] = x[i] + h

    return x, y


### PROBLEM PARAMETERS ###

CONNECTING_ROD_LENGHT = 15  # [cm]
CRANK_LENGHT = 5  # [cm]
R_C = CONNECTING_ROD_LENGHT / CRANK_LENGHT  # +/- = 3-3.5

COMPRESSION_RATIO = 10  # Gasoline : in [8,13]  /  Diesel : in [16,25]µ
HEAT_CAPACITY_RATIO = 1.3

HEAT_CAPACITY_gasoline = 2800  # [kJ/kg]
HEAT_CAPACITY_diesel = 2800  # [kJ/kg]
Q_tot = HEAT_CAPACITY_gasoline * 0.004 / 1000  # [J]

START_ANGLE = np.pi / 12  # in interval [15°,40°]
COMBUSTION_DURATION = np.pi / 6  # between 40° and 70°

ENGINE_DISPLACEMENT = 6e-4  # [m^3]


### Function considered : PRESSURE DIFFERENTIAL EVOLUTION ###

def dP(T, P_T):
    V_T = ENGINE_DISPLACEMENT * (0.5 * (1 - np.cos(T) + R_C - np.sqrt(R_C * R_C - np.sin(T) * np.sin(T))) +
                                 1 / (COMPRESSION_RATIO - 1))
    dV_dT = ENGINE_DISPLACEMENT / 2 * (np.sin(T) + 1 / (2 * np.sqrt(R_C * R_C - np.sin(T) * np.sin(T))) *
                                       2 * np.sin(T) * np.cos(T))
    dQ_dT = Q_tot / 2 * np.pi / COMBUSTION_DURATION * np.sin(np.pi * (T - START_ANGLE) / COMBUSTION_DURATION)
    dP_dT = - HEAT_CAPACITY_RATIO * P_T / V_T * dV_dT + (HEAT_CAPACITY_RATIO - 1) * 1 / V_T * dQ_dT
    return dP_dT


### EVOLUTION OF PRESSURE ###

step = 0.02
angles = np.linspace(-2 * np.pi, 2 * np.pi, 400)
P_init = 200000

Ts, Ps = RKF45('dP', P_init, angles, step)

plt.plot(Ts, Ps, '-r')
plt.title('Pressure in an engine cylinder', fontsize=18)
plt.xlabel('Angles   [RAD]', fontsize=14)
plt.ylabel('Pressure   [Pa]', fontsize=14)
#plt.savefig('Pressure.png', bbox_inches='tight')
plt.show()


