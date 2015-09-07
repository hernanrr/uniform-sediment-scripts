#!/usr/bin/env python
""" This script compares slope measurements by different methods

"""
from __future__ import division
import os
import re
import pdb
import numpy as np
from numpy.polynomial import polynomial as P
try:
    import cPickle as pickle
except:
    import pickle


home = os.path.expanduser("~")
flumepath = (home + '/Documents/Experiments/Data/1-input/flume')
sourcepath = (home + '/Documents/Experiments/Data/1-input/profiles')
outputpath = (home + '/Documents/Experiments/Data/3-output/profiles')
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

numbers = re.compile(r'(\d+)')
def NumericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


def load_pickle(fpickle):
    """Load the profiles from the pickle"""
    with open(fpickle, 'rb') as infile:
        pkl = pickle.load(infile)
    return pkl 

def compute_node_slope(y, x):
    """Computes the slope between two points
    Parameters:
    -----
    y: List, array or tuple of two or more values
    x : List, array or tupe of same dimensions as eta

    returns: S: List, array or tuple of dimension N-1 of eta and x,
    giving the slope(s)

    """
    # Create a container array for the slope values
    S = np.full_like(x, 0., dtype = float)
    # Compute the slope between the first and second nodes
    S[0] = - np.diff(y[0:2]) / np.diff(x[0:2])
    # compute the rest of the slopes using central differences
    for i, value in enumerate(y, start=1):
        try:
            # Compute the slopes using central differences
            S[i] = (y[i-1] - y[i+1]) / (x[i+1] - x[i-1])
        except IndexError:
            # Fail gracefully at the last node.
            pass
    # Compute the slopes between the delta front and toe.
    S[-1] =  (y[-2] - y[-1]) / (x[-1] - x[-2])
    return S

def compute_fit_slope(y, x):
    """Computes the slope of a linear regression through several points
    Parameters:
    -----------
    y: List, array or tuple of two or more values
    x : List, array or tupe of same dimensions as eta

    returns: S: List, array or tuple of dimension N-1 of eta and x,
    giving the slope(s)

    """
    _, m = P.polyfit(x, y, 1)
    return -m
    
def compute_slopes(x, eta, xi):
    """Computes the slopes"""
    # Compute the average 
    Sl_node = np.mean(compute_node_slope(eta, x))
    Sf_node = np.mean(compute_node_slope(xi, x))
    Sl_fit = compute_fit_slope(eta, x)
    Sf_fit = compute_fit_slope(xi, x)
    
    return Sl_node, Sl_fit, Sf_node, Sf_fit

def error(x1, x2):
    """Computes the relative error of two variables"""
    return x2/x1 - 1
    
    
def main():
    """Main routine for slope comparison"""
    print 'Script started'
    # Set how many nodes downstream to ignore. I suggest 2.
    ds_lim = 2
    # Add nodes upstream to ignore. Move all this to input file
    us_lim = 0
    # Load the profiles
    runs = ['equilibrium']#, 'aggradation']

    hdr = 'Gs, Q, date, sequence, Sl_avg, Sl_fit, Sf_avg, Sf_fit, % error Sl, % error Sf'
    for run in runs:
        # Create a dictionary to store all the results:
        table = []
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

        print 'Starting loop'
        for key in sorted(profiles, key=NumericalSort):
            x = profiles[key]['x']
            # Convert the measurements to meters for xi and eta
            xi = profiles[key]['wse'] / 100.
            eta = profiles[key]['bed'] / 100.
            #pdb.set_trace()
            Sl_avg, Sl_fit, Sf_avg, Sf_fit = compute_slopes(x[us_lim: -ds_lim],
                                                            eta[us_lim:
                                                                -ds_lim],
                                                            xi[us_lim:
                                                               -ds_lim])
            Gs = int(key.split('-')[0])
            Q = int(key.split('-')[1][0:2])
            date = int(key.split('-')[2])
            seq = int(key.split('-')[3])
            eSl = error(Sl_avg, Sl_fit)
            eSf = error(Sf_avg, Sf_fit)
            row = np.array([Gs, Q, date, seq, Sl_avg, Sl_fit, Sf_avg, Sf_fit,
                            eSl, eSf])
            table.append(row)
        output = np.vstack(table)
        fname = 'slopes.csv'
        print 'writing results to file'
        np.savetxt(fname, output, delimiter=',', header=hdr)

    print 'Script ended successfully'
    return

if __name__ == '__main__':
    main()
