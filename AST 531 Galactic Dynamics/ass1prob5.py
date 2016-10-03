from pylab import *
from numpy import linspace, sin, pi

def mkcurve(dist):
    angles = linspace(-180,180.0,400)
    return 8.0*sin(angles*pi/180)*(220/dist - 220/8.0)

ax = subplot(111)
ax.plot(linspace(-180,180,400),mkcurve(4),label='4 kpc')
ax.plot(linspace(-180,180,400),mkcurve(6),label='6 kpc')
ax.plot(linspace(-180,180,400),mkcurve(10),label='10 kpc')
ax.plot(linspace(-180,180,400),mkcurve(12),label='12 kpc')
legend(loc='upper left')
show()
