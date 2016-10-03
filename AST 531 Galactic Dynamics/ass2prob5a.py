import numpy as np
from numpy.linalg import norm
from matplotlib import pyplot

_nfw_dispersion = 200.0  # km/s
_nfw_radius = 5.0  # kpc
cell_size = 0.5  # kpc

cgs_G = 6.674e-8  # G in cgs (cm^3 g^-1 s^-2). Doesn't matter, factors out.

nfw_rho = _nfw_dispersion ** 2 / (4 * np.pi * cgs_G * _nfw_radius ** 2)


def nfw_potential(position):
    """
    Calculates Navarro-Frenk-White gravitational potential at a given position 
    in cell coordinates.
    Note: This is slow since the looping is not numpyized. But it's a homework
    problem so that's not particularly important
    @param position: 3-tuple cell index to calculate potential
    """
    position = np.array(position)  # Turn tuple into vector
    total = 0.0
    for cell in np.ndindex((61, 61, 61)):
        cell = np.array(cell)
        cell_dist = physical_dist(position, cell)
        # Only accumulate for cells that aren't the test cell
        if cell_dist > 0.0:
            total += nfw_density(cell) / cell_dist
    return -cgs_G * cell_size ** 3 * total


def nfw_density(position):
    """
    Calculates the Navarro-Frenk-White density at a given position in cell 
    coordinates
    @param position: 3-tuple cell index to calculate density
    """
    rad_mag = physical_dist(position, (30, 30, 30))
    if rad_mag > 0.0:
        return nfw_rho * _nfw_radius \
               / (rad_mag * (1 + rad_mag / _nfw_radius) ** 2)
    else:
        return 0.0
    

def analytic_nfw_potential(position):
    """
    Calculates the Navarro-Frenk-White potential analytically
    @param position: 3-tuple cell index to calculate potential 
    """
    rad_mag = physical_dist(np.array(position), (30, 30, 30))
    return -_nfw_dispersion ** 2 * np.log(1 + rad_mag / _nfw_radius) \
           * _nfw_radius / rad_mag


def physical_dist(cell1, cell2):
    """
    Finds the physical (kpc) distance between cell1 and cell2
    @param cell1: numpy array of first cell index
    @param cell2: numpy array of second cell index
    """
    return cell_size * norm(cell1 - cell2, 2)  # explicity use 2-norm


if __name__ == '__main__':

    radii = np.arange(30, 61)
    potentials = np.zeros(radii.shape, dtype=float)
    analytic_potentials = np.zeros(radii.shape, dtype=float)
    for i, rad in enumerate(radii):
        print('calculating radius %i...' % rad,)
        potentials[i] = nfw_potential((rad, 30, 30))
        analytic_potentials[i] = analytic_nfw_potential((rad, 30, 30))
        print('%2f (numeric) %2f (analytic)' %
              (potentials[i], analytic_potentials[i]))

# Plot everything up
# fig = pyplot.figure()
# ax = pyplot.subplot(111)
# ax.plot((radii-30)*cellSize,potentials, label='Numeric $\Phi$')
# ax.plot((radii-30)*cellSize,analPotentials, label='Analytic $\Phi$')
# ax.plot((radii-30)*cellSize,potentials-8650, label='Numeric $\Phi-8650$')
# ax.legend()
# pyplot.show()
