import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#differential eqn for the number of bacteria at time t
#T in kelvin
def logistic(N, t, T, Tmin, Tmax, b, c, Nmax):
    sqrtr = b*(T-Tmin)*(1-np.exp(c*(T-Tmax)))
    dNdt = ((sqrtr)**2)*N*(1 - (N/Nmax)) 
    return dNdt

#domain
t = np.linspace(0,3,100)

#parameters
T = 300 #room temp
Tmin = 276
Tmax = 322
b = 0.1
c = 0.1
Nmax = 10
N0 = 1 #initial num of bacteria

for i in range(5):
    y1 = odeint(logistic, N0, t, args=(T, Tmin, Tmax, b, c, Nmax))
    plt.plot(t, y1, label='T=' + str(T))
    T+=4.5

plt.axhline(y=Nmax, color='r', ls='--')
plt.xlabel('time')
plt.ylabel('num of bacteria')
plt.tight_layout()
plt.legend()
#plt.savefig('rays.pdf')
#plt.tight_layout()
plt.show()
plt.clf()


