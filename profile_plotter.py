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

# See http://matplotlib.org/users/pgf.html#font-specification
import matplotlib as mpl
# Set backend to PGF. This is good for embedding images in LaTeX
# Needs to be done before loading pyplot
mpl.use('pgf')
# Configure fonts and all else

from matplotlib import rc
from matplotlib import rcParams
#rc('font', size = 10, **{'family':'sans-serif','sans-serif':['Helvetica']})
## for palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['palatino']})
#rc('text', usetex=True)
rcParams['pdf.fonttype'] = 42 # Makes text be editable instead of being images
rcParams['text.color'] = 'black'

pgf_with_custom_preamble = {
    'font.size': 18, # set font size
    'font.family': 'serif', # use sans-serif/main font for text elements
    'text.usetex': True,    # use inline math for ticks
    'pgf.rcfonts': False,   # don't setup fonts from rc parameters
    'pgf.preamble': [
        r'\usepackage{siunitx}',
        # ... to force siunitx to actually use your fonts
        r'\sisetup{detect-all}',   
        r'\usepackage{metalogo}',
        # unicode math setup
        r'\usepackage[math-style=TeX,vargreek-shape=unicode]{unicode-math}',  
        r'\setmainfont{Helvetica Neue LT Pro 45 Light}', # serif font via preamble
    ]
}

mpl.rcParams.update(pgf_with_custom_preamble)

# Load pyplot after all configurations are set
import matplotlib.pyplot as plt


home = os.path.expanduser("~")
sourcepath = (home + '/Documents/Experiments/Data/1-input/profiles')
outputpath = (home + '/Documents/Experiments/Data/3-output/profiles')
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

def my_plotter(ax, xdata, ydata, param_dict):
    """
    A helper function to make a graph

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    xdata : array
       The x data

    ydata : array
       The y data

    param_dict : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    out = ax.plot(xdata, ydata, **param_dict)
    return out

def plot_profile(x, wse, bed, h, key):
    """Plots a profile to PDF, given coordinate, water surface elevation and
    bed elevation
    """
    # linear fits to the data
    b_wse_fit, m_wse_fit = P.polyfit(x, wse, 1)
    wse_fit = x * m_wse_fit + b_wse_fit
    b_bed_fit, m_bed_fit = P.polyfit(x, bed, 1) 
    bed_fit = x * m_bed_fit + b_bed_fit
    b_depth_fit, m_depth_fit = P.polyfit(x, h, 1)
    depth_fit = x * m_depth_fit + b_depth_fit    
    wse_equation = r'$\xi = {:.4f}x + {:.3f}$'.format(m_wse_fit, b_wse_fit)
    bed_equation = r'$\eta = {:.4f}x + {:.3f}$'.format(m_bed_fit, b_bed_fit)
    depth_equation = r'$h = {:.4f}x + {:.3f}$'.format(m_depth_fit, b_depth_fit)
    
    # Plot identification information
    run_name = key.split('-')
    run_type = run_name[-1]
    if run_type=='ag':
        run_type = 'aggradation'
    else:
        run_type = 'equilibrium'
        
    run_number = run_name[-2]
    run_flowrate = run_name[1][0:2]
    run_date = datetime.datetime.strptime(run_name[2],'%Y%m%d').date()
    feed_rate_units = r'\,\si{\g \per \minute}'
    flow_rate_units = r'\,\si{\l \per \s}'    
    feed_rate = r'G_s = {}{}'.format(run_name[0], feed_rate_units)
    flow_rate = r'Q = {}{}'.format(run_flowrate, flow_rate_units)
    # plot_title = '{} {} {} run {} of {}'.format(feed_rate, flow_rate,
    #                                             run_type, run_number,
    #                                             run_date)

    plot_title = '{}, {}'.format(feed_rate, flow_rate)

    plot_file_name = '{}.pdf'.format(key)

    # Create the figure
    fig = plt.figure(figsize=(18,8), tight_layout=True)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    # Specify shortcuts for some sane defaults:
    ls = '-' # linestyle
    #ls = 'None'
    ms = 6
    mec = 'None' # markeredgecolor
    mew = 0.1 # markeredgewidth
    m = u'^' #marker
    c = 'blue'  
    mfc = '#66c2a5'
    lw = 2

    # Create parameter dictionaries to use in the loops
    wse_params = {'ls': 'None', 'marker' : ur'v','ms' : ms, 'mfc' : 'blue',
                  'mec' : 'blue', 'mew' : 0.1, 'label': r'Water surface $(\xi)$'}
    bed_params = {'ls' : 'None', 'marker' : ur'^', 'ms' : ms, 'mfc' : 'black',
                  'mec' : 'black', 'mew' : 0.1, 'label':'Bed $(\eta)$'}
    depth_params = {'ls': 'None', 'marker' : ur'v','ms' : ms, 'mfc' : 'blue',
                    'mec' : 'blue', 'mew' : 0.1, 'label': r'Water depth $(h)$'}
    wse_fit_params = {'color': 'blue', 'ls': '-', 'lw': lw, 'label': wse_equation}
    bed_fit_params = {'color': 'black','ls': '-', 'lw': lw, 'label': bed_equation}
    depth_fit_params = {'color': 'blue', 'ls': '-', 'lw': lw, 'label': depth_equation}

    # Aggregate the plot parameters into lists, to facilitate the loop
    data = [wse_fit, bed_fit, wse, bed]
    params = [wse_fit_params, bed_fit_params, wse_params, bed_params]

    # Plot water surface elevation and bed elevation in a loop
    for ydata, param_dict in zip(data, params):
        my_plotter(ax1, x, ydata, param_dict)
    # Plot water depth
    data = [h, depth_fit]
    params = [depth_params, depth_fit_params]
    for ydata, param_dict in zip(data, params):
        my_plotter(ax2, x, ydata, param_dict)

    # ax1.legend(fontsize=12, loc='upper left', numpoints=3, handlelength=2,
    #           frameon=False)

    fig.suptitle(plot_title, fontsize = 20 )

    # Set the axes labels
    elev_units = r'\,\si{\per \meter}'
    xlabel = r'Downstream coordinate/(x{})'.format(elev_units)
    zlabel = r'Elevation/(z{})'.format(elev_units)
    hlabel = r'Depth/(h{})'.format(elev_units)

    # Format the figure
    for ax in fig.axes:
        ax.set_xlim(0,9)
        ax.set_ylim(0, 0.50)
        ax.xaxis.set_visible(True)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel(xlabel, fontsize = 20)
        if ax==ax1:
            ylabel = zlabel
        else:
            ylabel = hlabel 
        ax.set_ylabel(ylabel, fontsize = 20) 
        ax.tick_params(axis=u'both', labelcolor = 'black')   
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        ax.legend(fontsize = 18, loc='upper left', numpoints=3, handlelength=2,
                  frameon=False)
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
        for spine in spines_to_keep:
            ax.spines[spine].set_color('black')

    fig.savefig(plot_file_name, dpi=100, format='pdf', transparent=True,
                bbox_inches='tight', pad_inches=0.1, frameon=False) 
    print '{} writen to disk'.format(plot_file_name)
    plt.close('all')
    
    

def main():
    """Main routine"""
    print 'Script started' 
    # Load the profiles
    runs = ['equilibrium', 'aggradation']
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
#            if key.split('-')[0]=='1500':
            x = profiles[key]['x']
            wse = profiles[key]['wse'] / 100.
            bed = profiles[key]['bed'] / 100.
            h = wse - bed
            plot_profile(x, wse, bed, h, key)
#            pdb.set_trace()
#            else:
#                pass
              
              

    print 'Script completed successfully'

if __name__ == '__main__':
    main()
