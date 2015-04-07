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
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib import rcParams
rc('font', size = 10, **{'family':'sans-serif','sans-serif':['Helvetica']})
## for palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['palatino']})
rc('text', usetex=True)
rcParams['pdf.fonttype'] = 42
rcParams['text.color'] = 'black'
rcParams['text.latex.preamble'] = [
    r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    r'\sisetup{detect-all}',   # ... to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

numbers = re.compile(r'(\d+)')
def NumericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

home = os.path.expanduser("~")
sourcepath = (home + '/Documents/Experiments/Data/input/sonars')
outputpath = (home + '/Documents/Experiments/Data/output/sonars')
eq_in_path = os.path.join(sourcepath, 'equilibrium')
ag_in_path = os.path.join(sourcepath, 'aggradation')
eq_out_path = os.path.join(outputpath, 'equilibrium', 'plots')
ag_out_path = os.path.join(outputpath, 'aggradation', 'plots')
gp_source = (home + '/Documents/Experiments/Data/processed/profiles')


spines_to_remove = ['top', 'right']
spines_to_keep = ['bottom', 'left']

almost_black = '#262626'



def load_sonars(fpickle):
    """Load the sonar data from the pickle"""
    with open(fpickle, 'rb') as infile:
        pkl = pickle.load(infile)
    return pkl 


def plot_timeseries(d, h, run, name):
    """Plots sonar data time series to PDF, given coordinate, water surface
    elevation and bed elevation

    """
    # Labeling information
    run_name = run.split('-')
    run_type = run_name[-1]
    if run_type=='ag':
        run_type = 'aggradation'
    else:
        run_type = 'equilibrium'
#    run_number = run_name[2]
    run_date = datetime.datetime.strptime(run_name[1],'%Y%m%d').date()
    feed_rate_units = r'\,\si{\g \per \minute}'
    gs = np.int(run_name[0])
    feed_rate = r'\num{}{}'.format(gs, feed_rate_units)
    probs = 'probability distributions'
    plot_title = '{} {} run: Bed elevation {} and {}'.format(feed_rate, run_type,
                                                             name, probs) 
    
    plot_file_name = '{}_{}.pdf'.format(run, name)
    

    # Create a grid
    rows = 6
    columns = 2
    rwidth = [4, 1] # Width Ratio
    rheight = [1, 1] # Height Ratio
    gs = gridspec.GridSpec(rows, columns,
                           width_ratios=rwidth)
                           #height_ratios = rheight)  
    # Create figure
    fig = plt.figure(figsize=(12,8), tight_layout=True)

    #create subfigures
    # ax1 = fig.add_subplot(6,1,1,) # 6 rows, one column, first plot
    # ax2 = fig.add_subplot(6,1,2, sharex=ax1, sharey=ax1) 
    # ax3 = fig.add_subplot(6,1,3, sharex=ax1, sharey=ax1) 
    # ax4 = fig.add_subplot(6,1,4, sharex=ax1, sharey=ax1) 
    # ax5 = fig.add_subplot(6,1,5, sharex=ax1, sharey=ax1) 
    # ax6 = fig.add_subplot(6,1,6, sharex=ax1, sharey=ax1) 

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[2])
    ax3 = fig.add_subplot(gs[4])
    ax4 = fig.add_subplot(gs[6])
    ax5 = fig.add_subplot(gs[8])
    ax6 = fig.add_subplot(gs[10])

    ax7 = fig.add_subplot(gs[1])
    ax8 = fig.add_subplot(gs[3])
    ax9 = fig.add_subplot(gs[5])
    ax10 = fig.add_subplot(gs[7])
    ax11 = fig.add_subplot(gs[9])
    ax12 = fig.add_subplot(gs[11])

    ls = '-' # linestyle
    #ls = 'None'
    ms = 3
    mec = 'None' # markeredgecolor
    mew = 0.1 # markeredgewidth
    m = u'o' #marker
    c = 'gainsboro'  
    mfc = '#66c2a5'
    # Plot fluctuations around mean
    ax1.plot(d['t'], d['probe_6'], ms = ms, marker = m, ls = ls, color = c, 
             mfc=mfc, mec = mec, mew=mew)
    ax2.plot(d['t'], d['probe_5'], ms = ms, marker = m, ls = ls, color = c, 
             mfc=mfc,mec = mec, mew=mew)
    ax3.plot(d['t'], d['probe_4'], ms = ms, marker = m, ls = ls, color = c, 
             mfc=mfc,mec = mec, mew=mew)
    ax4.plot(d['t'], d['probe_2'], ms = ms, marker = m, ls = ls, color = c,
             mfc=mfc, mec = mec, mew=mew)
    ax5.plot(d['t'], d['probe_1'], ms = ms, marker = m, ls = ls, color = c,
             mfc=mfc, mec = mec, mew=mew)
    ax6.plot(d['t'], d['probe_3'], ms = ms, marker = m, ls = ls, color = c,
             mfc=mfc, mec = mec, mew=mew)

    ax7.plot(h['probe_6']['pe'], h['probe_6']['y'], ms = ms, marker = m, ls =
             ls, mfc=mfc, mec = mec, mew=mew, label='PDF')
    ax8.plot(h['probe_5']['pe'], h['probe_5']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc,mec = mec, mew=mew)
    ax9.plot(h['probe_4']['pe'], h['probe_4']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc,mec = mec, mew=mew)
    ax10.plot(h['probe_3']['pe'], h['probe_3']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc, mec = mec, mew=mew)
    ax11.plot(h['probe_2']['pe'], h['probe_2']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc, mec = mec, mew=mew)
    ax12.plot(h['probe_1']['pe'], h['probe_1']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc, mec = mec, mew=mew)
    
    # Plot the cummulative probability distribution function
    ax7.plot(h['probe_6']['Pe'], h['probe_6']['y'], ms = ms, marker = m, ls =
             ls,  mec = mec, mew=mew, label='CDF')
    ax8.plot(h['probe_5']['Pe'], h['probe_5']['y'], ms = ms, marker = m, ls =
             ls,  mec = mec, mew=mew)
    ax9.plot(h['probe_4']['Pe'], h['probe_4']['y'], ms = ms, marker = m, ls =
             ls,  mec = mec, mew=mew)
    ax10.plot(h['probe_3']['Pe'], h['probe_3']['y'], ms = ms, marker = m, ls =
             ls,   mec = mec, mew=mew)
    ax11.plot(h['probe_2']['Pe'], h['probe_2']['y'], ms = ms, marker = m, ls =
            ls,   mec = mec, mew=mew)
    ax12.plot(h['probe_1']['Pe'], h['probe_1']['y'], ms = ms, marker = m, ls =
              ls,   mec = mec, mew=mew)

    # Plot fake zero x-axis in all plots
    
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.plot([0,1], [0,0], ls='-', c = 'gainsboro', marker='None')

    spines_to_remove = ['top', 'right']
    spines_to_keep = ['bottom', 'left']

    # Plot whichever is smallest: one hour of data or full data.
    for ax in fig.axes:
        # try:
        #     ax.set_xbound(lower=d['t'][0], upper=d['t'][723])            
        # except IndexError:
        #     ax.set_xbound(lower=d['t'][0], upper=d['t'][-1])            
        # if name=='profiles':
        #     ax.set_ybound(lower=50, upper=300)            
        # else:
        ax.set_ybound(lower=-40, upper=40)
        ax.yaxis.set_ticks([-40, -20, 0, 20, 40])         
            
        ax.xaxis.set_visible(False)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.label.set_color(almost_black)
        ax.yaxis.label.set_color(almost_black)
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
        for spine in spines_to_keep:
            ax.spines[spine].set_color('black')

    # Format PDF and CDF axes
    for ax in [ax7, ax8, ax9, ax10, ax11, ax12]:
        ax.plot([0,1], [0,0], ls='-', c = 'gainsboro', marker='None')
        ax.set_xbound(lower=0.0, upper=1.0)
        ax.set_ybound(lower=-40, upper=40)
        ax.yaxis.set_ticks([-40, -20, 0, 20, 40])         
        ax.xaxis.set_ticks([0, 0.50,  1])                     
        ax.xaxis.set_visible(True)#False)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.label.set_color(almost_black)
        ax.yaxis.label.set_color(almost_black)
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
        for spine in spines_to_keep:
            ax.spines[spine].set_color('black')
        if ax==ax7:
            label = ax.set_ylabel('x = 3.68 m', fontsize = 12, color =
                                  'black', rotation = 'horizontal')
        elif ax==ax8:
            label = ax.set_ylabel('x = 4.68 m', fontsize = 12, color =
                                  'black', rotation = 'horizontal')
        elif ax==ax9:
            label = ax.set_ylabel('x = 5.68 m', fontsize = 12, color =
                                  'black', rotation = 'horizontal')
        elif ax==ax10:
            label = ax.set_ylabel('x = 6.68 m', fontsize = 12, color =
                                  'black', rotation = 'horizontal')
        elif ax==ax11:
            label = ax.set_ylabel('x = 7.68 m', fontsize = 12, color =
                                  'black', rotation = 'horizontal')
        elif ax==ax12:
            label = ax.set_ylabel('x = 8.68 m', fontsize = 12, color =
                                  'black', rotation = 'horizontal')
        ax.yaxis.set_label_coords(1.40, 0.425)


    # add labels and stuff
    ax6.xaxis.set_visible(True)
    ax12.xaxis.set_visible(True)
    fig.subplots_adjust(hspace=0)
    

    xlabel1 = r'$t / hours$'
    xlabel2 = r'$p$'
    
    ax6.set_xlabel(xlabel1, fontsize = 16, color = 'black',
                  rotation = 'horizontal')
    ax12.set_xlabel(xlabel2, fontsize = 16, color = 'black',
                  rotation = 'horizontal')

    ylabel1 = r'$ \left( \eta - \bar{\eta} \right) $ /mm '
    fig.text(0.00, 0.5, ylabel1, ha='center', va='center', rotation='vertical',
             fontsize = 16)

    ax7.legend(fontsize=10, frameon=False, bbox_to_anchor=(1.05, 1.25))
    fig.suptitle(plot_title, y= 1.05, fontsize=16, fontweight='bold')

#    elev_units = r'\,\si{\per \meter}'
#    xlabel = r'Downstream coordinate/(x{})'.format(elev_units)
#    ylabel = r'Elevation/(z{})'.format(elev_units)

    fig.savefig(plot_file_name, dpi=300, format='pdf',
                transparent=True, 
                bbox_inches='tight', 
                pad_inches=0.1,
                frameon=False) 
    print '{} written to disk'.format(plot_file_name)
    plt.close('all')
    

def plot_sigma(gp, d):
    """Plot sigma vs qbstar and stuff"""
    # Location of the sonars
    # x = [3.68, 4.68, 5.68, 6.68, 7.68, 8.68]
    # Create array for the loads
    taub_star = np.array([], dtype=float)
    tau_star_bs = np.array([], dtype=float)
    taub = np.array([], dtype=float)
    qb_star = np.array([], dtype=float)
    # Create empty arrays to store the sigmas, per probe
    s1 = np.array([], dtype=float)
    s2 = np.array([], dtype=float)
    s3 = np.array([], dtype=float)
    s4 = np.array([], dtype=float)
    s5 = np.array([], dtype=float)
    s6 = np.array([], dtype=float)
    # Create the x-axis variables
    # Store the sigmas-per-probe in arrays, by incremental feedrate
    for feedrate in sorted(gp, key= NumericalSort):
        # X-axis vectors
        tau_star_bs = np.append(tau_star_bs, np.mean(gp[feedrate]['Mean']['taub_star_s']) )
        taub = np.append(taub, np.mean(gp[feedrate]['Mean']['taub']) )
        tau_star_b = np.append(taub_star, np.mean(gp[feedrate]['Mean']['taub_star']) )
        qb_star = np.append(qb_star, np.mean(gp[feedrate]['Meta']['qb_star']) )
        # Y-axis vectors
        s1 = np.append(s1, d[feedrate]['probe_1']['sigma'])
        s2 = np.append(s2, d[feedrate]['probe_2']['sigma'])
        s3 = np.append(s3, d[feedrate]['probe_3']['sigma'])
        s4 = np.append(s4, d[feedrate]['probe_4']['sigma'])
        s5 = np.append(s5, d[feedrate]['probe_5']['sigma'])
        s6 = np.append(s6, d[feedrate]['probe_6']['sigma'])


    # Create figure
    fig = plt.figure(tight_layout=True)
    
    # Create a plotting axis
    ax1 = fig.add_subplot(111)
    
    # Prepare some formatting
    ls = '-' # linestyle
    #ls = 'None'
    ms = 4
    mec = 'None' # markeredgecolor
    mew = 0.1 # markeredgewidth
    m = u'o' #marker
    c = 'gainsboro'  
    mfc = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f',
           '#e5c494', '#b3b3b3']


    # Create labels
    l6 = r'\SI{3.68}{\m}'
    l5 = r'\SI{4.68}{\m}'
    l4 = r'\SI{5.68}{\m}'
    l2 = r'\SI{6.68}{\m}'
    l1 = r'\SI{7.68}{\m}'
    l3 = r'\SI{8.68}{\m}'

    # Plot
    x = tau_star_bs
    ax1.plot(x, s6, ms=ms, mec=mec, mew=mew, label=l6, ls=ls, c=c, marker=u'>',
             mfc=mfc[7])
    ax1.plot(x, s5, ms=ms, mec=mec, mew=mew, label=l5, ls=ls, c=c, marker=u'<',
             mfc=mfc[4])
    ax1.plot(x, s4, ms=ms, mec=mec, mew=mew, label=l4, ls=ls, c=c, marker=u'^',
             mfc=mfc[3])
    ax1.plot(x, s2, ms=ms, mec=mec, mew=mew, label=l2, ls=ls, c=c, marker=u'v',
             mfc=mfc[1])
    ax1.plot(x, s1, ms=ms, mec=mec, mew=mew, label=l1, ls=ls, c=c, marker=u'o',
             mfc=mfc[0])
    ax1.plot(x, s3, ms=ms, mec=mec, mew=mew, label=l3, ls=ls, c=c, marker=u's',
             mfc=mfc[2])
    # Format the figure/canvas
    spines_to_remove = ['top', 'right']
    spines_to_keep = ['bottom', 'left']
    for ax in fig.axes:
        #ax.set_xbound(lower=d['t'][0], upper=d['t'][723])            
        #ax.set_xbound(lower=d['t'][0], upper=d['t'][-1])            
        #ax.set_ybound(lower=-40, upper=40)
        #ax.yaxis.set_ticks([-40, -20, 0, 20, 40])         
            
        ax.xaxis.set_visible(False)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.label.set_color(almost_black)
        ax.yaxis.label.set_color(almost_black)
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
        for spine in spines_to_keep:
            ax.spines[spine].set_color('black')

    # add labels and stuff
    ax1.xaxis.set_visible(True)
    fig.subplots_adjust(hspace=0)

    # Add a legend
    ax1.legend(fontsize=10, loc='upper right', numpoints = 1, frameon=False,
               title='Probe location') 

    xlabel1 = r'$\tau^{*}_{\text{bs}}$'
    ylabel1 = r'$\sigma_e / \text{mm}$ '
    fig.text(0.5, 0.0, xlabel1, ha='center', va='center', rotation='horizontal',
             fontsize = 16)
    fig.text(0.00, 0.5, ylabel1, ha='center', va='center', rotation='vertical',
             fontsize = 16)

    plot_title = 'Variability of bed elevation fluctuations'
    fig.suptitle(plot_title, y= 1.02, fontsize=16, fontweight='bold')
    
    plot_file_name = 'sigma_plot.pdf'
    fig.savefig(plot_file_name, dpi=100, format='pdf',
                transparent=True, 
                bbox_inches='tight', 
                pad_inches=0.1,
                frameon=False) 
    print '{} written to disk'.format(plot_file_name)

    return


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
        # Get the pickles
        # first sonar
        f = run + '_sonars.pickle'
        sonars = load_sonars(f)
        # then general parameters
        f1 = os.path.join(gp_source, run, run + '_global_stats_summary.pickle')
        gp = load_sonars(f1)
        # Create a storage dictionary
        d = {}
        # Choose the output path
        if run=='equilibrium':
            os.chdir(eq_out_path)            
        else:
            os.chdir(ag_out_path)
        for run in sonars:
            # Give sensible names to the data
            t = sonars[run]['ContinuousSec'] / (60. * 60.)
            probe_1 = sonars[run]['Probe_1']
            probe_2 = sonars[run]['Probe_2']
            probe_3 = sonars[run]['Probe_3']
            probe_4 = sonars[run]['Probe_4']
            probe_5 = sonars[run]['Probe_5']
            probe_6 = sonars[run]['Probe_6']

            feedrate = run.split('-')[0]
            data = {feedrate:{}}
            data[feedrate] = {'t':t, 'probe_1':probe_1, 'probe_2':probe_2, 
                    'probe_3':probe_3, 'probe_4':probe_4, 'probe_5':probe_5,
                    'probe_6':probe_6}  
            

            fluctuations = {feedrate:{}}
            for key, value in data[feedrate].items():
                # No need to do anythig for the time information
                if key == 't':
                    fluctuations[feedrate].setdefault('t', t)
                else:
                    # Find the trend. Will be average for equilibrium. 
                    mu = np.mean(value)
                    # Detrend to get fluctuations around the mean. 
                    fluctuations[feedrate][key] = value - mu
            # Compute the histogram of the fluctuations
            histograms = {feedrate:{}}
            # A container to collect all values computed by this plotter
            d.setdefault(feedrate,{})
            for key, value in fluctuations[feedrate].items():
                if key == 't':
                    d[feedrate].setdefault('t', t)
                else:
                    bins = np.linspace(-50, 50, num=21)
                    hist, bin_edges = np.histogram(value, bins=bins,
                                                   density=False )
                    n = hist.sum()
                    pe = hist / (n + 1.)
                    Pe = 1 - pe.cumsum()
                    value[value<-50], value[value>50] = 0., 0. 
                    variance = np.sum( value ** 2. ) / ( n - 1. )
                    sigma = np.sqrt(variance)
                    y = ( bin_edges[0:-1] + bin_edges[1:] ) / 2.
                    histograms[feedrate][key] = {'hist':hist, 
                                                 'bin_edges':bin_edges, 
                                                 'pe':pe,
                                                 'Pe':Pe,
                                                 'variance':variance, 
                                                 'sigma':sigma, 'n':n, 'y':y}
                    d[feedrate].setdefault(key, {})
                    d[feedrate][key] = histograms[feedrate][key]

            plot_timeseries(fluctuations[feedrate], histograms[feedrate], run,
                            'fluctuations')    
            # Exits the loop for all feedrates
        plot_sigma(gp, d)



    print 'Script completed successfully'
    return

if __name__ == '__main__':
    main()
