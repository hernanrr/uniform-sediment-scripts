#!/usr/bin/env python
"""This script removes the sidewall effects from flume data.

"""

from __future__ import division
import numpy as np
import newton_raphson as nr
import pdb

g = np.float(9.81) # Gravity in m/s2
nu = np.float(1e-6) # Kinematic viscosity of water

rho = np.float(1.0e3) # Density of water, in kg/m3
rho_s = np.float(2.65e3) # Density of sediment, in kg/m3
R = (rho_s - rho) / rho # Submerged relative density of sediment

D = np.float(1.11e-3) # Characteristic diameter of the sediment in m


def fChezy(H, U, B, S):
    """Computes the Chezy friction coefficient"""
    return ( ( S * g ) / U ** 2 ) * ( (B * H) / ( B + 2 * H  ) )

def fRe(B, H, U, nu=1e-6):
    """Computes Reynolds number"""
    return ( U / nu ) * ( (B * H) / (B + 2 * H ) )

def fChezy_wall(xRef):
    """Computes the Chezy friction coefficient for the wall region"""
    # Do a first estimate of fw0 with the Blasius equation, as modified by
    # Chiew and Parker in 1994
    fw0 = 0.301 * xRef ** 0.2
    Rew0 = fw0 / xRef
    convergence = False
    while not convergence:
        Rew1 = Rew0
        fw0 = nr.newton_raphson(lambda fw: fNikuradse(fw, Rew0), 0.01)
        Rew0 = fw0 / xRef
        convergence = nr.good_enough(Rew1, Rew0)
    fw = fw0 / 8    
    return fw, Rew0


def fNikuradse(fw, Rew):
    """Nikuradse equation modified by for solving in a Newton-Raphson scheme.

    """
    nikuradse =  ( fw * ( 0.86 * np.log(4 * Rew * np.sqrt(fw) ) - 0.8 ) ** 2
                 - 1)
    return nikuradse


def remove_wall_effects(x, H, U, E, lim, B0=1.0):
    """Remove wall effects according to Vanoni and Brooks"""
    # Define some containers:
    # Area, bed-region
    Ab = np.full_like(x, 0., dtype=float)
    # Area, wall-region
    Aw = np.full_like(x, 0., dtype=float)
    # Chezy friction coefficient, total
    Cf = np.full_like(x, 0., dtype=float)
    # Chezy friction coefficient, bed-region
    Cfb = np.full_like(x, 0., dtype=float)
    # Chezy friction coefficient, wall-region
    Cfw = np.full_like(x, 0., dtype=float)
    # Reynolds number, total
    Re = np.full_like(x, 0., dtype=float)
    # Reynolds number, bed region
    Reb = np.full_like(x, 0., dtype=float)
    # Reynolds number, wall region
    Rew = np.full_like(x, 0., dtype=float)
    # Ratio of Chezy friction coefficient to Reynolds number
    Ref1 = np.full_like(x, 0., dtype=float)
    # Shear stress, bed-region
    taub = np.full_like(x, 0., dtype=float)
    # Shear stress, wall-region
    tauw = np.full_like(x, 0., dtype=float)
    # Sidewall-corrected Shields number, (bed-region)
    taub_star = np.full_like(x, 0., dtype=float)
    # Sidewall-corrected shear velocity, (bed-region)
    ub_star = np.full_like(x, 0., dtype=float)

    # Perform the computations Compute total energy slope. The -1
    # accounts for python's indexing quirks
    S = ( E[0] - E[lim-1] ) / ( x[lim-1] - x[0] )
    # Compute total friction coefficient
    Cf[:lim] = fChezy(H[:lim], U[:lim], B0, S)
    # Compute total Reynolds number
    Re[:lim] = fRe(B0, H[:lim], U[:lim], nu)
    # Compute the ratio of the friction coefficient to the
    # Reynolds No.
    Ref1[:lim] = Cf[:lim] / Re[:lim]
    # Compute wall-region friction coefficient, node per node
    # Get the wall-region Reynold's number while we are at it. 
    for i, value in enumerate(Cf[:lim]):
       Cfw[i], Rew[i] = fChezy_wall(Ref1[i] * 8)
    # Find the wall-region area
    Aw[:lim] = Rew[:lim] * nu * (2 * H[:lim] ) / U[:lim] 
    # Find the bed-region area by substrating wall-region area from
    # total area
    Ab[:lim] = B0 * H[:lim] - Aw[:lim]
    # Compute the bed-region friction factor
    Cfb[:lim] = Cf[:lim] + ( ( 2 * H[:lim] ) / B0 ) * ( Cf[:lim] -
                                                        Cfw[:lim] )
    # Compute the bed-region Reynolds number
    Reb[:lim] = Cfb[:lim] / Ref1[:lim]
    # Compute shear stresses for bed region
    taub[:lim] = rho * Cfb[:lim] * U[:lim] ** 2
    # Compute shear stresses for wall region
    tauw[:lim] = rho * Cfw[:lim] * U[:lim] ** 2
    # Compute the sidewall-corrected Shields number for the bed region
    taub_star[:lim] = Cfb[:lim] * U[:lim] ** 2 / ( R * g * D)
    # Compute the sidewall-corrected shear velocity
    ub_star[:lim] = np.sqrt( taub[:lim] / rho )
   
    # Collect the results
    r = Cf, Cfb, Cfw, Re, Reb, Rew, Ab, Aw, taub_star, taub, tauw, ub_star, S
    return r



# if __name__ == '__main__':
#     main()
