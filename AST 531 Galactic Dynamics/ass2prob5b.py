from __future__ import division
import numpy as np
from matplotlib import pyplot

nfw_dispersion = 200.0  # km/s
nfw_radius = 5.0  # kpc
cell_size = 0.5  # kpc

# I split this part into a separate file for my convenience
# Potentials, calculated from other file
potentials = np.array([-28789.53248398, -28012.72715114, -27189.71325385,
                       -25979.10634141, -24777.00779663, -23643.8818753,
                       -22589.70049077, -21611.80990653, -20704.16744864,
                       -19860.12727408, -19073.34885662, -18338.07476767,
                       -17649.17997044, -17002.1377608,  -16392.9584649,
                       -15818.12321584, -15274.52148289, -14759.39531239,
                       -14270.29083278, -13805.01660086, -13361.60802181,
                       -12938.29700704, -12533.48608081, -12145.72623491,
                       -11773.69792794, -11416.19471723, -11072.10909245,
                       -10740.42014816, -10420.18279014, -10110.51820834,
                       -9810.60063246]) - 8650


def nfw_radial_force(radius):
    radius = int(radius - 0.5)
    return -(potentials[radius + 1] - potentials[radius]) / cell_size


def nfw_v_circ(radius):
    return np.sqrt(-nfw_radial_force(radius) * radius * cell_size)


def analytic_nfw_v_circ(radius):
    radius = (radius-0.5) * cell_size  # Translate to physical units
    return nfw_dispersion * np.sqrt(np.log(1 + radius / nfw_radius)
                                    / (radius / nfw_radius)
                                    - 1 / (1 + radius / nfw_radius))

radii = np.linspace(0, 29, 30) + 0.5
velocities = np.zeros(radii.shape, dtype=float)
analytic_velocities = np.zeros(radii.shape, dtype=float)
for i, rad in enumerate(radii):
    velocities[i] = nfw_v_circ(rad)
    analytic_velocities[i] = analytic_nfw_v_circ(rad)

fig = pyplot.figure()
ax = pyplot.subplot(111)
ax.plot(radii * cell_size, velocities, label='Numeric $V_c$')
ax.plot(radii * cell_size, analytic_velocities, label='Analytic $V_c$')
ax.legend()
pyplot.show()
