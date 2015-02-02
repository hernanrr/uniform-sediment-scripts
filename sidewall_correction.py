#!/usr/bin/env python
""" This script removes the effects of the wall and form drag

"""
from __future__ import division
import os
import pdb
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle
import itertools


home = os.path.expanduser("~")
flumepath = (home + '/Documents/Experiments/Data/input/flume')
sourcepath = (home + '/Documents/Experiments/Data/input/profiles')
outputpath = (home + '/Documents/Experiments/Data/processed/profiles')
eq_in_path = os.path.join(sourcepath, 'equilibrium')
ag_in_path = os.path.join(sourcepath, 'aggradation')
eq_out_path = os.path.join(outputpath, 'equilibrium')
ag_out_path = os.path.join(outputpath, 'aggradation')

# Some constants. Should probably be read from a file as a dictionary...
g = np.float(9.81) # Gravity in m/s2
nu = np.float(1e-6) # Kinematic viscosity of water
rho = np.float(1.0e3) # Density of water, in kg/m3
B0 = np.float(19e-2) # Channel width in m
Q = np.float(30e-3) # water discharge in m3/s
xi_d = np.float(32e-2) # initial downstream tailgate elevation in m
D = np.float(1.11e-3) # Characteristic diameter of the sediment in m
rho_s = np.float(2.65e3) # Density of sediment, in kg/m3
R = (rho_s - rho) / rho # Submerged relative density of sediment 


def load_pickle(fpickle):
    """Load the profiles from the pickle"""
    with open(fpickle, 'rb') as infile:
        pkl = pickle.load(infile)
    return pkl 


def iterative_improve(good_enough, improve):
    """ """
#    pdb.set_trace()
    def iterate(guess):
        next_try = improve(guess)
        if good_enough(guess, next_try):
            return next_try
        else:
            return iterate(next_try)
    return lambda x: iterate(x)

def fixed_point(f, first_guess, epsilon=1.0e-8):
    """ """
    def good_enough(v1, v2):
        return abs(v1-v2)/v2 < epsilon
    return iterative_improve(good_enough, f)(first_guess)
        

def approx_deriv(f, dx=1.0e-8):
    """
    
    """
    return lambda x: float(f(x + dx) - f(x)) / dx

def newton_transform(f):
    """"""
    return lambda x: x - float(f(x)) / ( approx_deriv(f) (x) )
    
    
def newton_raphson(f, guess):
    """
    newton_raphson(f, guess)

    Newton Raphson procedure
    
    Parameters
    ----------
    f : function
        
    Returns
    -------
    float
        Zero of the specified function.
    """
    return fixed_point(newton_transform(f), guess)


def fU(Q, B, H):
    """Computes the average velocity per section"""
    return Q / (B0 * H)


def fFr(U, H):
    """Computes the Froude Number, given velocity and depth"""
    return U / np.sqrt(g * H)


def fE(xi, U):
    """Computes the specific energy for the cross-section"""
    return xi + U ** 2/(2*g)


def compute_friction_slope(E, x):
    """Computes the friction slope"""
    # Create a container for the friction slope
    Sf = np.full_like(x, 0., dtype = float)
    # Backward difference for the first node
    Sf[0] = (E[0] - E[1]) / (x[1] - x[0])
    # Central difference for the central nodes
    for i, value in enumerate(x, start=1):
        try:
            # Compute the slopes using central differences
            Sf[i] = (E[i-1] - E[i+1]) / (x[i+1] - x[i-1])
        except IndexError:
            # Fail gracefully at the extreme nodes (first and last)
            pass
    # Forward difference for the last node
    Sf[-1] = (E[-2] - E[-1]) / (x[-1] - x[-2])
    return Sf
    

def compute_bed_slope(eta, x):
    """Computes the slope between two points
    Parameters:
    -----
    eta: List, array or tuple of two or more values
    x : List, array or tupe of same dimensions as eta

    returns: S: List, array or tuple of dimension N-1 of eta and x,
    giving the slope(s)

    """
    # Create a container array for the slope values
    Sl = np.full_like(x, 0., dtype = float)
    # Compute the slope between the first and second nodes
    Sl[0] = - np.diff(eta[0:2]) / np.diff(x[0:2])
    # compute the rest of the slopes using central differences
    for i, value in enumerate(eta, start=1):
        try:
            # Compute the slopes using central differences
            Sl[i] = (eta[i-1] - eta[i+1]) / (x[i+1] - x[i-1])
        except IndexError:
            # Fail gracefully at the last node.
            pass
    # Compute the slopes between the delta front and toe.
    Sl[-1] =  (eta[-2] - eta[-1]) / (x[-1] - x[-2])
    return Sl

def fChezy(H, U, B, S):
    """Computes the Chezy friction coefficient"""
    return ( ( S * g ) / U ** 2 ) * ( (B * H) / ( B + 2 * H  ) ) 
    
def fRe(B, H, U, nu=1e-6):
    """Computes Reynolds number"""
    return ( U / nu ) * ( (B * H) / (B + 2 * H ) )

def fNikuradse(fw, Rew):
    """Nikuradse equation modified by for solving in a Newton-Raphson scheme.

    """
    nikuradse =  ( fw * ( 0.86 * np.log(4 * Rew * np.sqrt(fw) ) - 0.8 ) ** 2
                 - 1)
    return nikuradse


def good_enough(v1, v2, epsilon=1e-8):
    return abs( (v1 - v2) / v2 ) < epsilon



def fChezy_wall(xRef):
    """Computes the Chezy friction coefficient for the wall region"""
    # Do a first estimate of fw0 with the Blasius equation, as modified by
    # Chiew and Parker in 1994
    fw0 = 0.301 * xRef ** 0.2
    Rew0 = fw0 / xRef
    convergence = False
    while not convergence:
        Rew1 = Rew0
        fw0 = newton_raphson(lambda fw: fNikuradse(fw, Rew0), 0.01)
        Rew0 = fw0 / xRef
        convergence = good_enough(Rew1, Rew0)
    fw = fw0 / 8    
    return fw, Rew0
    
    

def remove_wall_effects(x, H, U, E, lim, B0=B0):
    """Remove wall effects according to Vanoni and Brookes"""
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
    # Compute total Chezy friction coefficient
    Cf[:lim] = fChezy(H[:lim], U[:lim], B0, S)
    # Compute total Reynolds number
    Re[:lim] = fRe(B0, H[:lim], U[:lim], nu)
    # Compute the ration of the Chezy friction coefficient to the
    # Reynolds No.
    Ref1[:lim] = Cf[:lim] / Re[:lim]
    # Compute wall-region Chezy friction coefficient, node per node
    # Get the wall-region Reynold's number while we are at it. 
    for i, value in enumerate(Cf[:lim]):
       Cfw[i], Rew[i] = fChezy_wall(Ref1[i] * 8)
    # Find the wall-region area
    Aw[:lim] = Rew[:lim] * nu * (2 * H[:lim] ) / U[:lim] 
    # Find the bed-region area by substrating wall-region area from
    # total area
    Ab[:lim] = B0 * H[:lim] - Aw[:lim]
    # Compute the bed-region Chezy friction factor
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
    
def compute_skin_friction(swc, lim):
    """Remove form drag effects from the side-wall-corrected results"""
    # Create a dictionary to store the result
    d = {}
    # Specify roughness height
    nk = np.float(3.5)
    # Specify coefficient to resistance relationship
    alpha_r = 8.1
    # Specify roughness height
    ks = nk * D
    # Define some vectors
    d['Rhb_s'] = np.full_like(swc['taub_star'], 0., dtype=float)
    d['ub_star_s'] = np.full_like(swc['taub_star'], 0., dtype=float)
    d['Cfbs'] = np.full_like(swc['taub_star'], 0., dtype=float)
    d['taub_star_s'] = np.full_like(swc['taub_star'], 0., dtype=float)
    d['phi'] = np.full_like(swc['taub_star'], 0., dtype=float)
    # Compute bed-region hydraulic radius
    d['Rhb'] = swc['Ab'] / B0
    for i, Rhb in enumerate(d['Rhb'][:lim]):
        convergence = False
        while not convergence:
            Rhb1 = Rhb
            ub_star_s = newton_raphson(lambda ub_star: fManning(ub_star,
                                                                swc['U'][i],
                                                                Rhb, alpha_r,
                                                                ks), 0.05)
            Rhb = ( ub_star_s ** 2 ) / (swc['S'] * g )
            convergence = good_enough(Rhb1, Rhb)
        # After convergence, store the converged values
        d['Rhb_s'][i] = Rhb
        d['ub_star_s'][i] = ub_star_s
        
    # We now choose the correct values:
    # First, specify the comparison condition
    condition = [ d['Rhb_s'][:lim] < d['Rhb'][:lim] ]
    # Then, specify the functions based on the conditions
    choice_Rhbs = [ d['Rhb_s'][:lim], d['Rhb'][:lim] ] 
    choice_ub_star_s = [ d['ub_star_s'][:lim], swc['ub_star'][:lim] ]
    choice_taub_star_s = [ d['ub_star_s'][:lim] ** 2 / ( R * g * D ) , \
                           swc['taub_star'][:lim] ]

    # Finally, apply the conditions and functions
    d['Rhb_s'][:lim] = np.where( condition, *choice_Rhbs )
    d['ub_star_s'][:lim] = np.where(  condition, *choice_ub_star_s )
    d['taub_star_s'][:lim] = np.where( condition, *choice_taub_star_s )
    d['phi'][:lim] =  d['taub_star_s'][:lim] / swc['taub_star'][:lim]
    d['Cfbs'][:lim] = ( d['ub_star_s'][:lim] / swc['U'][:lim] ) ** 2
    # We need Cfbs to plot the skin stresses. 
    return d, ks


def fManning(u_star_b, U, Rbs, alpha_r, ks):
    """Computes the Manning resistence relation, somehow"""
    manning =    U - u_star_b *  ( alpha_r * (Rbs / ks) ** (1/6) )
    return manning


def write_pickle(d, run, suffix):
    """Pickle the dictionary that is passed to the function"""
    # Dump the profile in the corresponding input folder
    wd = os.path.join(outputpath, run)
    os.chdir(wd)
    pickle_hdr = '{}_{}.pickle'.format(run, suffix)
    print 'Saving results in file {}/{}'.format(wd, pickle_hdr)
    with open(pickle_hdr, 'wb') as pickle_outfile:
        pickle.dump(d, pickle_outfile, -1) 
    return


def compute_statistics(d, lim):
    """Computes statistics of the variables stored in the dictionary"""
    stats={}
    for key in d:
        if key=='S':
            pass
        else:
            stats[key]={}
            # Count the number of non-zero elements in the array
            stats[key]['N'] = np.count_nonzero( d[key][:lim] )
            # Compute the arithmetic mean of the values
            stats[key]['Mean'] = np.mean( d[key][:lim] )
            # Compute the median of the values
            stats[key]['Median'] = np.median( d[key][:lim] )
            # Compute the Standard Deviation
            stats[key]['Std'] = np.std( d[key][:lim] )
            # Compute the variance
            stats[key]['Var'] = np.var( d[key][:lim] )
            # Compute the maximum
            stats[key]['Max'] = np.amax( d[key][:lim] )
            # Compute the minimum
            stats[key]['Min'] = np.amin( d[key][:lim] )
            # Compute the range (max - min)
            stats[key]['PtP'] = np.ptp( d[key][:lim] )
            # Compute the histograms
            stats[key]['Histogram'] = np.histogram( d[key][:lim] )
            # Consider adding the density argument to the histogram
    return stats
    

def summarize_statistics(stats):
    """Summarizes statistics per feedrate"""
    d={}
    for run, params in stats.iteritems():
        for param, stat_name in params.iteritems():
            for stat, value in stat_name.iteritems():
                d.setdefault(run.split('-')[0], {}
                ).setdefault(stat, {}
                ).setdefault(param, []
                ).append(value)
    return d

def add_feedrate_metadata(stats, ks):
    """Add feedrate data to the stats summary"""
    for run in stats:
        # Create the 'Meta' key and make the value be an empty dictionary
        stats.setdefault(run, {}).setdefault('Meta', {})
        # Put some LaTeX describing the feed rate
        stats[run]['Meta']['Gs'] = (r'{}{}'.format(run, '\,\si{\g \per \min}'))
        # Compute the volumetric unit feed rate
        stats[run]['Meta']['qf'] = np.int(run) / (1000 * B0 * 60 * rho_s )
        # Compute the Einstein Number for the run
        stats[run]['Meta']['qb_star'] = (stats[run]['Meta']['qf'] 
                                         / ( D * np.sqrt(R*g*D) ) )
        # Store the roughness height for each run
        stats[run]['Meta']['ks'] = ks
    return stats
    
    

def main():
    """Main routine for sidewall correction"""
    print 'Script started'
    # Specify the feedrates that this script considers for some things
    feed = ['500', '1000', '1500', '2000', '2500', '3000', '4000', '6000',
            '8000']
    # Set how many nodes downstream to ignore. I suggest 2.
    lim = 2
    # Load the profiles
    runs = ['equilibrium']
    for run in runs:
        # Create a dictionary to store all the results:
        d = {}
        # Create a dictionary to store statistics
        stats = {}
        # Choose source path
        if run=='equilibrium':
            os.chdir(eq_in_path)
        else:
            os.chdir(ag_in_path)
        # Get the pickle
        f = run + '_profiles.pickle'
        # Profiles is a dictionary of numpy structured arrays
        profiles = load_pickle(f)
        # Choose the output path
        if run=='equilibrium':
            os.chdir(eq_out_path)            
        else:
            os.chdir(ag_out_path)
        # Specify dictionary keys to collect side-wall corrected values
        swc_keys = ('Cf', 'Cfb', 'Cfw', 'Re', 'Reb', 'Rew', 'Ab', 'Aw',
                    'taub_star', 'taub', 'tauw', 'ub_star', 'S')

        for key in profiles:
            x = profiles[key]['x']
            # Convert the measurements to meters for xi and eta
            xi = profiles[key]['wse'] / 100.
            eta = profiles[key]['bed'] / 100.
            # Compute water depth
            H = xi - eta
            # Compute total area
            A = B0 * H
            # Compute the bed slope
            Sl = compute_bed_slope(eta, x)
            # Compute the velocities 
            U = fU(Q, B0, H)
            # Compute Froude number 
            Fr = fFr(U, H)
            # Compute the specific energies
            E = fE(xi, U)
            # Compute the friction slope
            Sf = compute_friction_slope(E, x)
            # Remove wall effects and collect results in a single variable.
            swc_values = remove_wall_effects(x, H, U, E, -lim, B0=B0)
            d[key] = {'x': x, 'xi': xi, 'eta': eta, 'H': H, 'A': A, 'Sl': Sl,
            'U': U, 'Fr': Fr, 'E': E, 'Sf': Sf}
            d[key].update( dict(itertools.izip(swc_keys, swc_values)) )
            skin_friction, ks = compute_skin_friction(d[key], -lim)
            d[key].update(skin_friction)
            # Compute statistics on values stored in dictionary
            statistics = compute_statistics(d[key], -lim)
            stats[key]={}
            stats[key].update(statistics)
            # Summarize statistics
            stats_summary = summarize_statistics(stats)
            # Compute the mass feed rate
            d[key]['Gs'] = np.int(key.split('-')[0])
            # Put it in a nice LaTeX string
            d[key]['Feed rate'] = (r'{:d}{}' .format(
                d[key]['Gs'],'\,\si{\g \per \minute}' ) )
            # Compute the volumetric unit feed rate
            d[key]['qf'] = d[key]['Gs'] / (1000 * B0 * 60 * rho_s )
            # Compute the Einstein Number for the run
            d[key]['qb_star'] = d[key]['qf'] / (D * np.sqrt( R * g * D ) )
            # Store the roughness height
            d[key]['ks'] = ks
        # Add metadata do the stats summary (i.e. Einstein's number)
        stats_summary = add_feedrate_metadata(stats_summary, ks)
        # Store all the input parameters
        d['input'] = {'g':g, 'nu':nu, 'rho':rho, 'B0':B0, 'Q':Q,
                      'xi_d':xi_d, 'D':D, 'rho_s':rho_s, 'R':R}
        # Once all the values are computed, pickle the results:
        write_pickle(d, run, 'global_parameters')
        # Pickle the statistics
        write_pickle(stats, run, 'global_statistics')
        # Pickle the statistics summary
        write_pickle(stats_summary, run, 'global_stats_summary')
    print 'Script completed successfully'
    return
if __name__ == '__main__':
    main()
