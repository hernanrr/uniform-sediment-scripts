#!/usr/bin/env python

"""This script reads the profiles stored in the "input" folder as pickles and
plots them in the output file

"""
from __future__ import division
import glob
#import sys
import os
import pdb
import re
import datetime
import numpy as np
from numpy.polynomial import polynomial as P
try:
    import cPickle as pickle
except:
    import pickle

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
rc('font', size = 10, **{'family':'sans-serif','sans-serif':['Helvetica']})
## for palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['palatino']})
rc('text', usetex=True)
rcParams['pdf.fonttype'] = 42
rcParams['text.color'] = 'dimgray'
rcParams['text.latex.preamble'] = [
    r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]


home = os.path.expanduser("~")
sourcepath = (home + '/Documents/Experiments/Data/input/profiles')
outputpath = (home + '/Documents/Experiments/Data/output/profiles')
eq_in_path = os.path.join(sourcepath, 'equilibrium')
ag_in_path = os.path.join(sourcepath, 'aggradation')
eq_out_path = os.path.join(outputpath, 'equilibrium', 'plots')
ag_out_path = os.path.join(outputpath, 'aggradation', 'plots')

spines_to_remove = ['top', 'right']
spines_to_keep = ['bottom', 'left']

almost_black = '#262626'



def load_profiles(fpickle):
    """Load the profiles from the pickle"""
    with open(fpickle, 'rb') as infile:
        pkl = pickle.load(infile)
    return pkl 

def plot_profile(x, wse, bed, key):
    """Plots a profile to PDF, given coordinate, water surface elevation and
    bed elevation"""
    # linear fits to the data
    b_wse_fit, m_wse_fit = P.polyfit(x, wse, 1)
    wse_fit = x * m_wse_fit + b_wse_fit
    b_bed_fit, m_bed_fit = P.polyfit(x, bed, 1) 
    bed_fit = x * m_bed_fit + b_bed_fit
    wse_equation = r'$\xi = {:.4f}x + {:.3f}$'.format(m_wse_fit, b_wse_fit)
    bed_equation = r'$\eta = {:.4f}x + {:.3f}$'.format(m_bed_fit, b_bed_fit)
    run_name = key.split('-')
    run_type = run_name[-1]
    if run_type=='ag':
        run_type = 'aggradation'
    else:
        run_type = 'equilibrium'
    run_number = run_name[2]
    run_date = datetime.datetime.strptime(run_name[1],'%Y%m%d').date()
    feed_rate_units = r'\,\si{\g \per \minute}'
    feed_rate = '{}{}'.format(run_name[0], feed_rate_units)
    plot_title = '{} {} run {} of {}'.format(feed_rate, run_type, run_number,
                                             run_date) 
    plot_file_name = '{}.pdf'.format(key)
    fig = plt.figure(figsize=(12,3), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.scatter(x, wse, s=8, facecolor = 'blue', edgecolor = 'blue', lw=0,
               label= r'Water surface $(\xi)$', marker= ur'v' )
    ax.plot(x, wse_fit, color='blue', linestyle='-', label=wse_equation)
    ax.scatter(x, bed, s=8,  facecolor = 'grey', edgecolor = 'grey', lw=0,
               label='Bed $(\eta)$')
    ax.plot(x, bed_fit, color='gray', linestyle='-', label=bed_equation)
    ax.legend(fontsize=12, loc='upper left', frameon=False)
    fig.suptitle(plot_title, fontsize=12)

    # Set the axes labels
    elev_units = r'\,\si{\per \meter}'
    xlabel = r'Downstream coordinate/(x{})'.format(elev_units)
    ylabel = r'Elevation/(z{})'.format(elev_units)

    # Format the figure
    for ax in fig.axes:
        ax.set_xlim(0,9)
        ax.set_ylim(0, 0.50)
        ax.xaxis.set_visible(True)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel(xlabel, fontsize = 12)
        ax.set_ylabel(ylabel, fontsize = 12) 
        ax.tick_params(axis=u'both', labelcolor = 'gray')   
        ax.xaxis.label.set_color('gray')
        ax.yaxis.label.set_color('gray')
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    for spine in spines_to_keep:
        ax.spines[spine].set_color('gray')
    fig.savefig(plot_file_name, dpi=100, format='pdf',
                transparent=True, bbox_inches='tight', pad_inches=0.1,
                frameon=False) 
    print '{} writen to disk'.format(plot_file_name)
    plt.close('all')
    
    

def main():
    """Main routine"""
    print 'Script started' 
    # Load the profiles
    runs = ['equilibrium']
    for run in runs:
        # Choose source path
        if run=='equilibrium':
            os.chdir(eq_in_path)            
        else:
            os.chdir(ag_in_path)
        # Get the pickle
        f = run + '_profiles.pickle'
        profiles = load_profiles(f)
        # Choose the output path
        if run=='equilibrium':
            os.chdir(eq_out_path)            
        else:
            os.chdir(ag_out_path)
        for key in profiles:
#            if key.split('-')[0]=='16000':
            x = profiles[key]['x']
            wse = profiles[key]['wse'] / 100.
            bed = profiles[key]['bed'] / 100.
            plot_profile(x, wse, bed, key)             
#            else:
#                pass
              
              

    print 'Script completed successfully'


if __name__ == '__main__':
    main()
