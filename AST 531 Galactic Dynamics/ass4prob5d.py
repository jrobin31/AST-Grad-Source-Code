from __future__ import division, print_function
import numpy as np
from matplotlib import pyplot

step_size = 1e-3
steps = int(10 / step_size)

phi0 = 0.0
vel_r0 = 0.0
vel_z0 = 0.0
rad0 = 0.31
height0 = 0.17
q = 0.9

# Angular momenta of orbits we will calculate
Lz = 0.2


# Equation for dPhi/dt in Axisymmetric potential
def phi_dot(ang_momentum, radius):
    return ang_momentum / radius ** 2


# Equation for dv_r/dt in Axisymmetric potential
def vel_r_dot(ang_momentum, radius, height):
    return ang_momentum ** 2 / radius ** 3 \
           - radius / (radius ** 2 + height ** 2)


# Equation for dv_z/dt in Axisymmetric potential
def vel_z_dot(q, radius, height):
    return -height / q / (radius**2 + height**2)

# Using leapfrog integration, find velocity at t=1/2
vel_r_half = vel_r0 + step_size / 2. * vel_r_dot(Lz, rad0, height0)
vel_z_half = vel_z0 + step_size / 2. * vel_z_dot(q, rad0, height0)

# Rows are timesteps, columns are values of L
phis = np.zeros((steps,))
radii = np.zeros((steps,))
heights = np.zeros((steps,))
vel_rs = np.zeros((steps,))
vel_zs = np.zeros((steps,))

# Set initial conditions for each column (L value)
phis[0] = phi0
radii[0] = rad0
heights[0] = height0
vel_zs[0] = vel_z_half
vel_rs[0] = vel_r_half

for step in range(1, steps):
    phis[step] = phis[step-1] + step_size * phi_dot(Lz, radii[step - 1])
    radii[step] = radii[step-1] + step_size * vel_rs[step - 1]
    heights[step] = heights[step-1] + step_size * vel_zs[step - 1] / q
    vel_rs[step] = vel_rs[step - 1] + step_size * vel_r_dot(
        Lz, radii[step], heights[step])
    vel_zs[step] = vel_zs[step - 1] + step_size * vel_z_dot(
        q, radii[step], heights[step])

    next_pct = 100 * (step + 1) // steps
    curr_pct = 100 * step // steps
    if next_pct - curr_pct > 0:
        print('{:d}%'.format(next_pct))
print('')

# find zero crossings
up_zc_mask = np.diff(np.sign(heights), 1) > 0

plot_title = r'$\~v_{r,\tau=0} = %0.2g$, $\~v_{z,\tau=0} = %0.2g$,' \
             r'$\~r_{\tau=0} = %0.2g$, $\~z_{\tau=0} = %0.2g$' % \
             (vel_r0, vel_z0, rad0, height0)

# Plot r vs. phi
fig = pyplot.figure()
ax = pyplot.subplot(111)
ax.plot(radii, heights, label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
# pyplot.xlim(0.0,0.5)
ax.set_xlabel('$\~r$')
ax.set_ylabel('$\~z$')
ax.legend(loc='lower right')
ax.set_title(plot_title)
pyplot.show()

# Plot y vs. x
ax = pyplot.subplot(111, aspect='equal', adjustable='datalim')
x = radii*np.cos(phis)
y = radii*np.sin(phis)
ax.plot(x, y, label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.legend(loc='upper right')
ax.set_title(plot_title)
pyplot.show()

# Plot upward zero crossings, r vs. v_r
ax = pyplot.subplot(111)
ax.plot(radii[1:][up_zc_mask], vel_rs[1:][up_zc_mask], 'o',
        label='$\~L_z = %0.2g$, $q = %0.2g$' % (Lz, q))
ax.set_xlabel('$\~r$')
ax.set_ylabel('$\~v_r$')
ax.legend(loc='lower right')
ax.set_title(plot_title)
pyplot.show()
