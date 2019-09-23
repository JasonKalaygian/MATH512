import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import scipy as sp

#datapoints (webplotdigitizer)
ecolipts = np.array([[293.1672377451225, 5.2035203740147615],
[296.06685114332237, 7.633597394683473],
[299.0943328371777, 7.96068049908594],
[301.1999640370418, 8.93653536857037],
[303.22667652318114, 9.804081795750374],
[305.093953228175, 10.331298762274853],
[307.19199224798456, 10.719518895537593],
[309.3683505988831, 11.169655255087257],
[311.1567085901521, 11.58856377930731],
[313.0159935266675, 11.49721786559843],
[315.03211691957284, 11.545168476469232],
[317.2028810325365, 11.562310819855545],
[318.0350239253569, 9.970170724154118],
[319.26755441894846, 9.368030928144009],
[321.21454901451506, 4.065412624996251]])

ecolitemp = ecolipts[:, 0]
ecolisqr = 1e1 * ecolipts[:, 1]

#sqrt(rate of growth) function
def sqrtr(T, Tmin, Tmax, b, c):
    return b*(T-Tmin)*(1-np.exp(c*(T-Tmax)))

def ecoli(T, b, c):
    return b*(T-276)*(1-np.exp(c*(T-322)))

model = sp.optimize.curve_fit(ecoli, ecolitemp, ecolisqr)
print(model)
bparam = model[0][0]
cparam = model[0][1]

#lsq = sp.optimize.least_squares(ecoli, 

#differential eqn for the number of bacteria at time t
#T in kelvin
def logistic(N, t, T, Tmin, Tmax, b, c, Nmax, N0, f):
    sqr = sqrtr(T, Tmin, Tmax, b, c)
    Nmin = (1-1e-6)*N0
    r = sqr**2/3600
    #print(r)
    dNdt = (r)*N*(1 - (N/Nmax))*(1 - (Nmin/N))**f
    return dNdt

#domaini
t = np.linspace(0,10, 10000)

#parameters
T = 300 #room temp
Tmin = 276
Tmax = 322
b = bparam
c = cparam
Nmax = 1e9
N0 = 1e2 #initial num of bacteria
f=0.72

for i in range(5):
    y1 = odeint(logistic, N0, t, args=(T, Tmin, Tmax, b, c, Nmax, N0, f))
    plt.semilogy(t, y1, label='T=' + str(T))
    T+=4.5

plt.axhline(y=Nmax, color='r', ls='--')
plt.xlabel('time [hrs]')
plt.ylabel('num of bacteria $N$')
plt.tight_layout()
plt.legend()
plt.show()
plt.clf()

Ts = np.linspace(Tmin-5, Tmax+5)
plt.plot(Ts, ecoli(Ts, bparam, cparam))
#plt.plot(ecolitemp, ecolisqr, 'o')
plt.vlines(314, 0, 130, linestyles='dashed', colors='r')
plt.xlabel('Temperature [K]')
plt.ylabel('$\sqrt{r}$')
plt.title('E. Coli')
plt.ylim(0, 130)
plt.savefig('growthcurve.pdf')
plt.show()
plt.clf()
