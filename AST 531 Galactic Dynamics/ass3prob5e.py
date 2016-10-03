import numpy as np
from matplotlib import pyplot
from progressbar import ProgressBar

stepSize = 1e-3
steps = 20000

phiStart = 0.0
rStart = 2 - 1e-6

# Angular momenta of orbits we will calculate
Lc = 2*np.sqrt(np.log(3)/2 - 1/3.)
angMomentaDivisors = np.array([4., 2., 1.5, 1.2, 1.1, 1.05]) 
angMomenta = Lc/angMomentaDivisors


# Equation for dPhi/dt in NFW
def phiDot(ang_momentum, radius):
    return ang_momentum / radius ** 2


# Equation for dr/dt in NFW
def rDot(ang_momentum, radius):
    return -np.sqrt(2 * (np.log(1 + radius) / radius - np.log(3) / 2 - ang_momentum ** 2 / 2 * (1 / radius ** 2 - 1 / 4.)))


## Rows are timesteps, columns are values of L
phis = np.zeros((steps,len(angMomenta)))
radii = np.zeros((steps,len(angMomenta)))

## Set initial conditions for each column (L value)
phis[0, :] = phiStart
radii[0, :] = rStart

## Add a progress bar in the terminal so we know how far along the calculation is
pBar = ProgressBar(maxval=len(phis.flat)).start()

for amStep, angMomentum in enumerate(angMomenta):
    for step in xrange(1,steps):
        radii[step,amStep] = radii[step-1,amStep] + stepSize*rDot(angMomentum, radii[step-1,amStep])
        phis[step,amStep] = phis[step-1,amStep] + stepSize*phiDot(angMomentum, radii[step-1,amStep])
        pBar.update(step*amStep)
pBar.finish()

## Plot r vs. phi
fig = pyplot.figure()
ax = pyplot.subplot(111)
for amStep,divisor in enumerate(angMomentaDivisors):
    ax.plot(phis[:,amStep], radii[:,amStep], label='$\~L = \~L_c/%0.3g$' % divisor)
ax.axvline(np.pi)
ax.axvline(np.pi/2)
ax.set_xlabel('$\phi$')
ax.set_ylabel('$\~r$')
ax.set_title('NFW Orbital Parameters, apocenter $\~r=2$')
ax.legend(loc='lower right')
pyplot.show()

## Plot y vs. x
ax = pyplot.subplot(111, aspect='equal')
for amStep,divisor in enumerate(angMomentaDivisors):
    x = radii[:,amStep]*np.cos(phis[:,amStep])
    y = radii[:,amStep]*np.sin(phis[:,amStep])
    ax.plot(x, y, label='$\~L = \~L_c/%0.3g$' % divisor)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('NFW Orbits, apocenter $\~r=2$')
ax.legend(loc='upper right')
pyplot.show()