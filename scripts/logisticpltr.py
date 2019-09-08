import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def logistic(N, t, T, Tmin, Tmax, b, c, Nmax):
    sqrtr = b*(T-Tmin)*(1-np.exp(c*(T-Tmax)))
    dNdt = ((sqrtr)**2)*N*(1 - (N/Nmax)) 
    return dNdt

t = np.linspace(0,2,100)

y1 = odeint(logistic, 1, t, args=(300, 276, 322, 0.1, 0.1, 10))
plt.plot(t, y1)

plt.xlabel('time')
plt.ylabel('num of bacteria')
plt.tight_layout()
#plt.savefig('rays.pdf')
#plt.tight_layout()
plt.show()
plt.clf()


