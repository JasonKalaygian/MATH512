import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#sqrt(rate of growth) function
def sqrtr(T, Tmin, Tmax, b, c):
    return b*(T-Tmin)*(1-np.exp(c*(T-Tmax)))

#differential eqn for the number of bacteria at time t
#T in kelvin
def logistic(N, t, T, Tmin, Tmax, b, c, Nmax, N0, f):
    sqr = sqrtr(T, Tmin, Tmax, b, c)
    Nmin = (1-1e-6)*N0
    dNdt = ((sqr)**2)*N*(1 - (N/Nmax))*(1-(Nmin/N))**f 
    return dNdt

#domain
t = np.linspace(0,5,10000)

#parameters
T = 300 #room temp
Tmin = 276
Tmax = 322
b = 0.1
c = 0.1
Nmax = 1e6
N0 = 1 #initial num of bacteria

'''
for i in range(5):
    y1 = odeint(logistic, N0, t, args=(T, Tmin, Tmax, b, c, Nmax, N0, 0.8))
    plt.semilogy(t, y1, label='T=' + str(T))
    #plt.ylim([0, 6.5])
    T+=4.5
'''

y1 = odeint(logistic, N0, t, args=(T, Tmin, Tmax, b, c, Nmax, N0, 0.8))
plt.semilogy(t, y1, label='T=' + str(T))

plt.axhline(y=Nmax, color='r', ls='--')
plt.xlabel('time $t$')
plt.ylabel('$N$')
plt.tight_layout()
#plt.legend()
#plt.show()
#plt.clf()
plt.savefig('sigmoid.pdf')

'''
Ts = np.linspace(Tmin-5, Tmax+5)
plt.plot(Ts, sqrtr(Ts, Tmin, Tmax, b, c))
plt.ylim(0, 3)
plt.xlabel('Temperature [K]')
plt.ylabel('sqrt(r)')
plt.show()
plt.clf()
'''
