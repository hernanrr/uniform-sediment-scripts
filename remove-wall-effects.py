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

import sidewall_correction as sw
import newton_raphson as nr
from const import g, nu


home = os.path.expanduser("~")
flumepath = (home + '/Documents/Experiments/Data/1-input/flume')
sourcepath = (home + '/Documents/Experiments/Data/1-input/profiles')
outputpath = (home + '/Documents/Experiments/Data/2-processed/profiles')
eq_in_path = os.path.join(sourcepath, 'equilibrium')
ag_in_path = os.path.join(sourcepath, 'aggradation')
eq_out_path = os.path.join(outputpath, 'equilibrium')
ag_out_path = os.path.join(outputpath, 'aggradation')

# Some constants. Should probably be read from a file as a dictionary...

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
    

        
def compute_Einstein_skin_friction(swc, ds_lim):
    """Remove form drag effects from the side-wall-corrected results"""
    # Create a dictionary to store the result
    d = {}
    # Specify roughness height
    nk = np.float(3.5) #4.5:6, 6.4:8, 9.1:10
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
    for i, Rhb in enumerate(d['Rhb'][:ds_lim]):
        convergence = False
        while not convergence:
            Rhb1 = Rhb
            ub_star_s = nr.newton_raphson(lambda ub_star:
                                       fManning(ub_star,
                                                swc['U'][i],
                                                Rhb,
                                                alpha_r,
                                                ks),
                                       0.05)
            Rhb = ( ub_star_s ** 2 ) / (swc['S'] * g )
            convergence = nr.good_enough(Rhb1, Rhb)
        # After convergence, store the converged values
        d['Rhb_s'][i] = Rhb
        d['ub_star_s'][i] = ub_star_s
        
    # We now choose the correct values:
    # First, specify the comparison condition
    condition = [ d['Rhb_s'][:ds_lim] < d['Rhb'][:ds_lim] ]
    # Then, specify the functions based on the conditions
    choice_Rhbs = [ d['Rhb_s'][:ds_lim], d['Rhb'][:ds_lim] ] 
    choice_ub_star_s = [ d['ub_star_s'][:ds_lim], swc['ub_star'][:ds_lim] ]
    choice_taub_star_s = [ d['ub_star_s'][:ds_lim] ** 2 / ( R * g * D ) , \
                           swc['taub_star'][:ds_lim] ]

    # Finally, apply the conditions and functions
    d['Rhb_s'][:ds_lim] = np.where( condition, *choice_Rhbs )
    d['ub_star_s'][:ds_lim] = np.where(  condition, *choice_ub_star_s )
    d['taub_star_s'][:ds_lim] = np.where( condition, *choice_taub_star_s )
    d['phi'][:ds_lim] =  d['taub_star_s'][:ds_lim] / swc['taub_star'][:ds_lim]
    d['Cfbs'][:ds_lim] = ( d['ub_star_s'][:ds_lim] / swc['U'][:ds_lim] ) ** 2
    # We need Cfbs to plot the skin stresses. 
    return d, ks

def Engelund_Hansen_skin_friction(swc, ds_lim):
    """Remove form drag effects from the side-wall-corrected results using the
    Engelund-Hansen decomposition

    Formulation
    -----------
    \begin{equation}
        \theta' = 0.4\theta^2
    \end{equation}
    
    Where:
    """
    
    # Create a dictionary to store the result
    d = {}
    # Define some vectors
    d['EH_tau_star_b_s'] = np.full_like(swc['taub_star'], 0., dtype=float)
    
    d['EH_tau_star_b_s'][:ds_lim] = 0.4 * swc['taub_star'][:ds_lim] ** 2 
    return d

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


def compute_statistics(d, ds_lim):
    """Computes statistics of the variables stored in the dictionary"""
    stats={}
    for key in d:
        if key=='S':
            pass
        else:
            stats[key]={}
            # Count the number of non-zero elements in the array
            stats[key]['N'] = np.count_nonzero( d[key][:ds_lim] )
            # Compute the arithmetic mean of the values
            stats[key]['Mean'] = np.mean( d[key][:ds_lim] )
            # Compute the median of the values
            stats[key]['Median'] = np.median( d[key][:ds_lim] )
            # Compute the Standard Deviation
            stats[key]['Std'] = np.std( d[key][:ds_lim] )
            # Compute the variance
            stats[key]['Var'] = np.var( d[key][:ds_lim] )
            # Compute the maximum
            stats[key]['Max'] = np.amax( d[key][:ds_lim] )
            # Compute the minimum
            stats[key]['Min'] = np.amin( d[key][:ds_lim] )
            # Compute the range (max - min)
            stats[key]['PtP'] = np.ptp( d[key][:ds_lim] )
            # Compute the histograms
            stats[key]['Histogram'] = np.histogram( d[key][:ds_lim] )
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
    # Set how many nodes downstream to ignore. I suggest 2.
    ds_lim = 2
    # Add nodes upstream to ignore. Move all this to input file
    us_lim = 0
    # Load the profiles
    runs = ['equilibrium', 'aggradation']
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
            x = profiles[key]['x'][us_lim:-ds_lim] 
            # Convert the measurements to meters for xi and eta
            xi = profiles[key]['wse'][us_lim:-ds_lim] / 100.
            eta = profiles[key]['bed'][us_lim:-ds_lim] / 100.
            # Compute water depth
            H = np.mean(xi - eta)
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

            swc_values = sw.remove_wall_effects(x, H, U, E, -ds_lim, B0=B0)
            # Create dictionary to store the computed global parameters
            d[key] = {'x': x, 'xi': xi, 'eta': eta, 'H': H, 'A': A, 'Sl': Sl,
            'U': U, 'Fr': Fr, 'E': E, 'Sf': Sf}
            # Add side-wall corrected values to the global parameter dictionary
            d[key].update( dict(itertools.izip(swc_keys, swc_values)) )

            
            # Remove form drag from the sidewall-corrected parameters. Using
            # Einstein decomposition
            Einstein_skin_friction, ks = compute_Einstein_skin_friction(d[key], -ds_lim)
            # Update the global-parameters dictionary with skin friction values.
            d[key].update(Einstein_skin_friction)
            # Remove form drag from the siewall-corrected parameters using
            # Engelund-Hansen decomposition.
            Engelund_Hansen_friction = Engelund_Hansen_skin_friction(d[key],
                                                                     ds_lim)
            # Update the global-parameters dictionary with skin friction values.
            d[key].update(Engelund_Hansen_friction)

            # Compute statistics on values stored in dictionary
            statistics = compute_statistics(d[key], -ds_lim)
            # Create a dictionary to store statistics and fill it. 
            stats[key]={}
            stats[key].update(statistics)
            # Summarize statistics across runs of same feedrate. 
            stats_summary = summarize_statistics(stats)
            # Compute the mass feed rate
            d[key]['Gs'] = np.int(key.split('-')[0])
            # Put it in a nice LaTeX string
            d[key]['Feed rate'] = (r'{:d}{}' .format(
                d[key]['Gs'],'\,\si{\g \per \minute}' ) )
            # Compute the volumetric unit feed rate, convert grams and minutes.
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
        #pdb.set_trace()
    print 'Script completed successfully'
    return
if __name__ == '__main__':
    main()
