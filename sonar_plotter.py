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

import matplotlib as mpl
mpl.use('pgf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
#from matplotlib import font_manager
from matplotlib import rcParams
##rc('font', size = 10, **{'family':'sans-serif','sans-serif':['Helvetica']})
## for palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['palatino']})
#rc('text', usetex=True)
rcParams['pdf.fonttype'] = 42 # Makes text be editable instead of being images
rcParams['text.color'] = 'black'

pgf_with_custom_preamble = {
    'font.size': 10, # set font size
    "font.family": "serif", # use sans-serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
         "\\usepackage{siunitx}",         # load additional packages
        r'\sisetup{detect-all}',   # ... to force siunitx to actually use your fonts
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",  # unicode math setup
         r"\setmathfont{xits-math.otf}",
         r"\setmainfont{HelveticaNeueLTPro-Roman}", # serif font via preamble
         ]
}

mpl.rcParams.update(pgf_with_custom_preamble)
# rcParams['text.latex.preamble'] = [
#     r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#     r'\sisetup{detect-all}',   # ... to force siunitx to actually use your fonts
#     r'\usepackage{helvet}',    # set the normal font here
#     r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#     r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
# ]

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

# Wavelengths, in meters
eq_lambda = {'500': 0.860, '1000': 0.860, '1500': 0.773, '2000': 0.620, '2500':
              0.520, '3000': 0.460, '4000': 0.378, '6000': 0.418, '8000':
              0.427, '10000': 0.317, '12000': 0, '16000': 0, '20000': 0 }

# Amplitudes, in cm
eq_amplitude = {'500': 4.3, '1000': 3.8, '1500': np.nan, '2000'
                 : np.nan, '2500': np.nan, '3000': 4, '4000': 2.75, '6000': 3.85
                , '8000': 3.825, '10000': np.nan, '12000': 4.633, '16000':
                 np.nan , '20000 ' : np.nan}

# H_nod, in meters
eq_h_nod = {'500': 0.1263, '1000': 0.1235, '1500': 0.1471, '2000': 0.1209,
            '2500': 0.1219, '3000': 0.1094, '4000': 0.1038, '6000': 0.1039,
            '8000': 0.0916, '10000': 0.0938, '12000': 0.0940, '16000': 0.0853,
            '20000': 0.0858}

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
    flow_rate_units = r'\,\si{\l \per \s}'
    gs = np.int(run_name[0])
    qw = np.int(run_name[2][0:2])
    feed_rate = r'\num{}{}'.format(gs, feed_rate_units)
    flow_rate = r'\num{}{}'.format(qw, flow_rate_units)
    probs = 'probability distributions'
    plot_title = '{}, {} {} run: Bed elevation {} and {}'.format(feed_rate,
                                                                 flow_rate,
                                                                 run_type, name,
                                                                 probs) 
    
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
    ax10.plot(h['probe_2']['pe'], h['probe_2']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc, mec = mec, mew=mew)
    ax11.plot(h['probe_1']['pe'], h['probe_1']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc, mec = mec, mew=mew)
    ax12.plot(h['probe_3']['pe'], h['probe_3']['y'], ms = ms, marker = m, ls =
             ls,  mfc=mfc, mec = mec, mew=mew)
    
    # Plot the cummulative probability distribution function
    ax7.plot(h['probe_6']['Pe'], h['probe_6']['y'], ms = ms, marker = m, ls =
             ls,  mec = mec, mew=mew, label='CDF')
    ax8.plot(h['probe_5']['Pe'], h['probe_5']['y'], ms = ms, marker = m, ls =
             ls,  mec = mec, mew=mew)
    ax9.plot(h['probe_4']['Pe'], h['probe_4']['y'], ms = ms, marker = m, ls =
             ls,  mec = mec, mew=mew)
    ax10.plot(h['probe_2']['Pe'], h['probe_2']['y'], ms = ms, marker = m, ls =
             ls,   mec = mec, mew=mew)
    ax11.plot(h['probe_1']['Pe'], h['probe_1']['y'], ms = ms, marker = m, ls =
            ls,   mec = mec, mew=mew)
    ax12.plot(h['probe_3']['Pe'], h['probe_3']['y'], ms = ms, marker = m, ls =
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


def plot_pdf_comparison(h, runtype):
    """Plots a comparison of the PDF and CDF of the bed elevation fluctuations
    by plotting, respectively, the PDF and CDF of a few feedrates in a single
    plot. 

    """
    # Figure size
    figsize = (4,4) # (width, height) in inches because of the 100 dpi thing
    # below in the savefig function.

    # Some formatting shortcuts
    ls = '-' # linestyle
    lw = 1
    ms = 3 # markersisze
    mec = 'None' # markeredgecolor
    mew = 0.1 # markeredgewidth
    m = u'o' #marker
    c = 'gainsboro' # light-gray color for fake axis 
    #mfc = '#66c2a5'
    
    # Comment out the undesired runs.
    runs_to_plot = [
        '500',
        '1000',
        '1500',
        #'2000',
        #'2500',
        '3000',
        #'4000',
        '6000',
        '8000',
        '10000',
        '12000',
        #'16000',
        '20000'
        ]

    # Units for the legend
    feed_rate_units = r'\,\si{\g \per \minute}'
    flow_rate_units = r'\,\si{\l \per \s}'

    # Create a figure.
    fig1 = plt.figure(figsize=figsize, tight_layout=True)
    fig2 = plt.figure(figsize=figsize, tight_layout=True)
    fig3 = plt.figure(figsize=figsize, tight_layout=True)
    fig4 = plt.figure(figsize=figsize, tight_layout=True)
    fig5 = plt.figure(figsize=figsize, tight_layout=True)
    fig6 = plt.figure(figsize=figsize, tight_layout=True)
    fig7 = plt.figure(figsize=figsize, tight_layout=True)
    fig8= plt.figure(figsize=figsize, tight_layout=True)
    fig9 = plt.figure(figsize=figsize, tight_layout=True)
    fig10 = plt.figure(figsize=figsize, tight_layout=True)
    fig11 = plt.figure(figsize=figsize, tight_layout=True)
    fig12 = plt.figure(figsize=figsize, tight_layout=True)

    

    # Create axis to plot the PDF and CDF
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ax3 = fig3.add_subplot(111)
    ax4 = fig4.add_subplot(111)
    ax5 = fig5.add_subplot(111)
    ax6 = fig6.add_subplot(111)
    ax7 = fig7.add_subplot(111)
    ax8 = fig8.add_subplot(111)
    ax9 = fig9.add_subplot(111)
    ax10 = fig10.add_subplot(111)
    ax11 = fig11.add_subplot(111)
    ax12 = fig12.add_subplot(111)




    # Populate the figure
    for flowrate in sorted(h, key=NumericalSort):
        for feedrate in sorted(h[flowrate], key=NumericalSort):
            if feedrate in runs_to_plot:
                print 'plotting {}'.format(feedrate)
                label = r'\num{}{}\,\num{}{}'.format(feedrate, feed_rate_units, flowrate,
                                           flow_rate_units)
                # Plot the cummulative probability distribution function
                ax1.plot(h[flowrate][feedrate]['probe_1']['Pe'], h[flowrate][feedrate]['probe_1']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax3.plot(h[flowrate][feedrate]['probe_2']['Pe'], h[flowrate][feedrate]['probe_2']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax5.plot(h[flowrate][feedrate]['probe_3']['Pe'], h[flowrate][feedrate]['probe_3']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax7.plot(h[flowrate][feedrate]['probe_4']['Pe'], h[flowrate][feedrate]['probe_4']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax9.plot(h[flowrate][feedrate]['probe_5']['Pe'], h[flowrate][feedrate]['probe_5']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax11.plot(h[flowrate][feedrate]['probe_6']['Pe'], h[flowrate][feedrate]['probe_6']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)




                # Plot the probability distribution function
                ax2.plot(h[flowrate][feedrate]['probe_1']['pe'], h[flowrate][feedrate]['probe_1']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax4.plot(h[flowrate][feedrate]['probe_2']['pe'], h[flowrate][feedrate]['probe_2']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax6.plot(h[flowrate][feedrate]['probe_3']['pe'], h[flowrate][feedrate]['probe_3']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax8.plot(h[flowrate][feedrate]['probe_4']['pe'], h[flowrate][feedrate]['probe_4']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax10.plot(h[flowrate][feedrate]['probe_5']['pe'], h[flowrate][feedrate]['probe_5']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)
                ax12.plot(h[flowrate][feedrate]['probe_6']['pe'], h[flowrate][feedrate]['probe_6']['y'],
                         ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                         label=label)         

            else:
                pass
        


    # Define axes labels
    xlabel_cdf = r'$p$'
    ylabel_cdf = r"$ \eta' = \left( \eta - \bar{\eta} \right) $ /mm "
    xlabel_pdf = r'$p$'
    ylabel_pdf = r"$ \eta' = \left( \eta - \bar{\eta} \right) $ /mm "

    for ax in [ax1, ax3, ax5, ax7, ax9, ax11]:
        ax.set_xlabel(xlabel_cdf, fontsize = 10, color = 'black',
                       rotation = 'horizontal')
        ax.set_ylabel(ylabel_cdf, fontsize = 10, color = 'black',
                       rotation = 'vertical')
    for ax in [ax2, ax4, ax6, ax8, ax10, ax12]:        
        ax.set_xlabel(xlabel_pdf, fontsize = 10, color = 'black',
                       rotation = 'horizontal')
        ax.set_ylabel(ylabel_pdf, fontsize = 10, color = 'black',
                       rotation = 'vertical')


    # Format figures
    figures = [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10,
               fig11, fig12]
    for fig in figures:
        for ax in fig.axes:
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
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
            
    for ax in axes:
        ax.plot([0,1], [0,0], ls='-', c = 'gainsboro', marker='None')
        ax.set_xbound(lower=0.0, upper=1.05)
        ax.set_ybound(lower=-40, upper=40)
        ax.yaxis.set_ticks([-40, -20, 0, 20, 40])         
        ax.xaxis.set_ticks([0, 0.50,  1])                     
        ax.xaxis.set_visible(True)#False)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
        for spine in spines_to_keep:
            ax.spines[spine].set_color('black')
        ax.legend(fontsize=10, loc='upper left', bbox_to_anchor = (1.05, 1), frameon=False)
        ax.tick_params(axis='both', colors='black')


    # Specify the titles
    #plot_title_cdf = 'CDF of bed elevation fluctuations'
    #plot_title_pdf = 'PDF of bed elevation fluctuations'
    #fig1.suptitle(plot_title_cdf, y= 1.05, fontsize=16)#, fontweight='bold')
    #fig2.suptitle(plot_title_pdf, y= 1.05, fontsize=16)#, fontweight='bold')

    # Save the file
    for fig in figures:
        if fig==fig1:
            probe = 'probe_1'
        elif fig==fig2:
            probe = 'probe_1'
        elif fig==fig3:
            probe = 'probe_2'
        elif fig==fig4:
            probe = 'probe_2' 
        elif fig==fig5:
            probe = 'probe_3'           
        elif fig==fig6:
            probe = 'probe_3'
        elif fig==fig7:
            probe = 'probe_4'            
        elif fig==fig8:
            probe = 'probe_4'
        elif fig==fig9:
            probe = 'probe_5'            
        elif fig==fig10:
            probe = 'probe_5'
        elif fig==fig11:
            probe = 'probe_6'             
        elif fig==fig12:
            probe = 'probe_6'            
        if fig in [fig1, fig3, fig5, fig7, fig9, fig11]:               
            plot_file_name = 'cdf_comparison_{}_{}.pdf'.format(probe, runtype)
        elif fig in [fig2, fig4, fig6, fig8, fig10, fig12]:
            plot_file_name = 'pdf_comparison_{}_{}.pdf'.format(probe, runtype)
            
        fig.savefig(plot_file_name, dpi=100, format='pdf',
                     transparent=True,
                     bbox_inches='tight',
                     pad_inches=0.1,
                     frameon=False)
    # Confirm that file has been saved
    print '{} written to disk'.format(plot_file_name)
        
    # Close the figure
    plt.close('all')

    ######################################
    # # Plot at all probes per feed rate #
    ######################################

    
    for feedrate in sorted(h, key=NumericalSort):
        if feedrate in runs_to_plot:
            fig13 = plt.figure(figsize=figsize, tight_layout=True)
            fig14 = plt.figure(figsize=figsize, tight_layout=True)
            ax13 = fig13.add_subplot(111)
            ax14 = fig14.add_subplot(111)        
            print 'plotting {}'.format(feedrate)

            ax13.plot(h[feedrate]['probe_6']['Pe'], h[feedrate]['probe_6']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 3.68 m')
            ax13.plot(h[feedrate]['probe_5']['Pe'], h[feedrate]['probe_5']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 4.68 m')
            ax13.plot(h[feedrate]['probe_4']['Pe'], h[feedrate]['probe_4']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 5.68 m')
            ax13.plot(h[feedrate]['probe_2']['Pe'], h[feedrate]['probe_2']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 6.68 m')
            ax13.plot(h[feedrate]['probe_1']['Pe'], h[feedrate]['probe_1']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 7.68 m')
            ax13.plot(h[feedrate]['probe_3']['Pe'], h[feedrate]['probe_3']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 8.68 m')

            
            ax14.plot(h[feedrate]['probe_6']['pe'], h[feedrate]['probe_6']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 3.68 m')
            ax14.plot(h[feedrate]['probe_5']['pe'], h[feedrate]['probe_5']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 4.68 m')
            ax14.plot(h[feedrate]['probe_4']['pe'], h[feedrate]['probe_4']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 5.68 m')
            ax14.plot(h[feedrate]['probe_2']['pe'], h[feedrate]['probe_2']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 6.68 m')
            ax14.plot(h[feedrate]['probe_1']['pe'], h[feedrate]['probe_1']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 7.68 m')
            ax14.plot(h[feedrate]['probe_3']['pe'], h[feedrate]['probe_3']['y'],
                     ms = ms, marker = m, ls = ls, lw = lw, mec = mec, mew=mew,
                      label='x = 8.68 m')
            # Define axes labels
            ######################################## is p the same for cdf and pdf?
            xlabel_cdf = r'$p$'
            ylabel_cdf = r"$ \eta' = \left( \eta - \bar{\eta} \right) $ /mm "
            xlabel_pdf = r'$p$'
            ylabel_pdf = r"$ \eta' = \left( \eta - \bar{\eta} \right) $ /mm "

            for ax in [ax13]:
                ax.set_xlabel(xlabel_cdf, fontsize = 12, color = 'black',
                               rotation = 'horizontal')
                ax.set_ylabel(ylabel_cdf, fontsize = 12, color = 'black',
                               rotation = 'vertical')
            for ax in [ax14]:
                ax.set_xlabel(xlabel_pdf, fontsize = 12, color = 'black',
                               rotation = 'horizontal')
                ax.set_ylabel(ylabel_pdf, fontsize = 12, color = 'black',
                               rotation = 'vertical')


            # Format figures
            figures = [fig13, fig14]
            for fig in figures:
                for ax in fig.axes:
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
            axes = [ax13, ax14]
            for ax in axes:
                ax.plot([0,1], [0,0], ls='-', c = 'gainsboro', marker='None')
                ax.set_xbound(lower=0.0, upper=1.05)
                ax.set_ybound(lower=-40, upper=40)
                ax.yaxis.set_ticks([-40, -20, 0, 20, 40])         
                ax.xaxis.set_ticks([0, 0.50,  1])                     
                ax.xaxis.set_visible(True)#False)
                ax.xaxis.set_ticks_position('none')
                ax.yaxis.set_ticks_position('none')
                ax.xaxis.label.set_color('black')
                ax.yaxis.label.set_color('black')
                for spine in spines_to_remove:
                    ax.spines[spine].set_visible(False)
                for spine in spines_to_keep:
                    ax.spines[spine].set_color('black')
                ax.legend(fontsize=10, frameon=False)
                ax.tick_params(axis='both', colors='black')


            # Specify the titles
            #plot_title_cdf = 'CDF of bed elevation fluctuations'
            #plot_title_pdf = 'PDF of bed elevation fluctuations'
            #fig13.suptitle(plot_title_cdf, y= 1.05, fontsize=16)#, fontweight='bold')
            #fig14.suptitle(plot_title_pdf, y= 1.05, fontsize=16)#, fontweight='bold')

            # Save the file
            figures2 = [fig13, fig14]
            for fig in figures2:           
                if fig==fig13:               
                    plot_file_name = 'x_cdf_plot_netherlands_{}_{}.pdf'.format(feedrate, runtype)
                else:
                    plot_file_name = 'x_pdf_plot_netherlands_{}_{}.pdf'.format(feedrate,
                                                                             runtype)
                fig.savefig(plot_file_name, dpi=100, format='pdf',
                            transparent=True,
                            bbox_inches='tight',
                            pad_inches=0.1,
                            frameon=False)
                # Confirm that file has been saved
                print '{} written to disk'.format(plot_file_name)

            # Close the figure
            plt.close('all')        

        else:
            pass
        
    return



def main():
    """Main routine"""
    print 'Script started'
    # Load the profiles
    runs = ['equilibrium']#, 'aggradation']
    for run_type in runs:
        # Choose source path
        if run_type=='equilibrium':
            os.chdir(eq_in_path)            
        else:
            os.chdir(ag_in_path)
        # Get the pickles
        # first sonar
        f = run_type + '_sonars.pickle'
        sonars = load_sonars(f)
        # then general parameters
        f1 = os.path.join(gp_source, run_type, run_type + '_global_stats_summary.pickle')
        gp = load_sonars(f1)
        # Create a storage dictionary
        d = {}
        probability = {}
        # Choose the output path
        if run_type=='equilibrium':
            os.chdir(eq_out_path)
            pickle_path = eq_in_path
        else:
            os.chdir(ag_out_path)
            pickle_path = ag_in_path
        # Aggregate the data by flowrate then feedrate instead of by
        # filename. This removes the date field stored in the `sonars`
        # dictionary. Start by creating a new dictionary
        data = {}  # Different for equilibrium and aggradation runs.
        # New dictionary as dictionary key to store fluctuation data
        fluctuations = {} 
        histograms = {}

        for run in sorted(sonars, key=NumericalSort):
            # Give sensible names to the data
            t = sonars[run]['RunSeconds'] / (60. * 60.)
            T = sonars[run]['ContinuousSec'] / (60. * 60.)
            probe_1 = sonars[run]['Probe_1']
            probe_2 = sonars[run]['Probe_2']
            probe_3 = sonars[run]['Probe_3']
            probe_4 = sonars[run]['Probe_4']
            probe_5 = sonars[run]['Probe_5']
            probe_6 = sonars[run]['Probe_6']
            
            # Get the feedrate information for the run we are working on
            feedrate = run.split('-')[0]
            flowrate = run.split('-')[2][:2]
            runtype = run.split('-')[-1]

            # Store non-detrended data
            data.setdefault(flowrate, {}).setdefault(feedrate, {})
            data[flowrate][feedrate] = {'t':t, 'T': T, 'probe_1':probe_1,
                                        'probe_2':probe_2, 'probe_3':probe_3,
                                        'probe_4':probe_4, 'probe_5':probe_5,
                                        'probe_6':probe_6}  
            # Store detrended data
            fluctuations.setdefault(flowrate, {}).setdefault(feedrate, {})

            # Detrending procedure, given instantanues values in data dict.
            for key, value in data[flowrate][feedrate].items():
                # Ensure that time starts at zero
                if key == 't':
                    fluctuations[flowrate][feedrate].setdefault('t', t - t[0])
                elif key == 'T':
                    fluctuations[flowrate][feedrate].setdefault('T', T - T[0])
                else:
                    # Find the trend. Assume very low slope for equilibrium
                    # This allows for the same procedure in aggradation
                    # index of np.nan values in value-array:
                    idx = np.isfinite(value)
                    b_bed_fit, m_bed_fit = P.polyfit((T-T[0])[idx], value[idx],
                                                     1)
                    bed_fit = (T-T[0]) * m_bed_fit + b_bed_fit
                    bed_equation = r'$\eta = {:.4f}x + {:.3f}$'.format(m_bed_fit, b_bed_fit)
                    #print bed_equation
                    mu = np.nanmean(value)
                    # Detrend to get fluctuations around the mean.
                    # Uncomment next line for fluctuation around zero
                    #fluctuations[feedrate][key] = value - mu
                    # Uncomment next time for fluctuation around a trend.
                    fluctuations[flowrate][feedrate][key] = value - bed_fit
            # Compute the histogram of the fluctuations
            histograms.setdefault(flowrate, {}).setdefault(feedrate, {})
            # A container to collect all values computed by this plotter
            d.setdefault(flowrate, {}).setdefault(feedrate, {})
            probability.setdefault(flowrate, {}).setdefault(feedrate, {})
            for key, value in fluctuations[flowrate][feedrate].items():
                if key == 't':
                    d[flowrate][feedrate].setdefault('t', t-t[0])
                elif key == 'T':
                    d[flowrate][feedrate].setdefault('T', T-T[0])                    
                else:
                    # Establish bounds for the bins and number of bins. This
                    # comes out to about 5 mm per bin.
                    idx = np.isfinite(value)
                    bins = np.linspace(-50, 50, num=21)
                    hist, bin_edges = np.histogram(value, bins=bins,
                                                   density=False )
                    n = np.nansum(hist)
                    pe = hist / (n + 1.)
                    Pe = 1 - pe.cumsum()
                    # Zero out values outside of the bin-bounds. This shouldn't
                    # be necessary after the sonar data is cleaned. We live
                    # with it for now. EDIT: Sonar data has been cleaned.
                    value[value[idx]<-50], value[value[idx]>50] = 0., 0. 
                    variance = np.nansum( value ** 2. ) / ( n - 1. )
                    sigma = np.sqrt(variance)
                    y = ( bin_edges[0:-1] + bin_edges[1:] ) / 2.
                    histograms[flowrate][feedrate][key] = {'hist':hist,
                                                           'bin_edges':bin_edges,
                                                           'pe':pe, 'Pe':Pe,
                                                           'variance':variance,
                                                           'sigma':sigma,
                                                           'n':n, 'y':y}
                    d[flowrate][feedrate].setdefault(key, {})
                    d[flowrate][feedrate][key] = histograms[flowrate][feedrate][key]

                    
            # plot_timeseries(fluctuations[flowrate][feedrate],
            #                 histograms[flowrate][feedrate], run,
            #                 'fluctuations')
            
        # Save the data for either the equilibrium or aggradational run
        pickle_names = ['histograms', 'fluctuations', 'data', 'd']
        pickle_contents = [histograms, fluctuations, data, d]

        for pickle_name, pickle_content in zip(pickle_names, pickle_contents):
            pickle_hdr = os.path.join(pickle_path, run_type + '_' +
                                      pickle_name) + '.pickle'  
            with open(pickle_hdr, 'wb') as pickle_outfile:
                pickle.dump(pickle_content, pickle_outfile, -1)

        for qw in sorted(histograms, key=NumericalSort):
            for qf in sorted(histograms[qw], key=NumericalSort):
                for probe in sorted(histograms[qw][qf], key=NumericalSort):
                    s = histograms[qw][qf][probe]['sigma']
                    #hist = histograms[qw][qf][probe]['hist']
                    pe = histograms[qw][qf][probe]['pe']
                    y = histograms[qw][qf][probe]['y']                    
                    print '{},{},{},{}'.format(qf, qw, probe, s)
                    #print '{},{}'.format(pe, y)            
        #     # Exits the loop for all feedrates
        # plot_sigma(gp, d)
        plot_pdf_comparison(d, runtype)


    print 'Script completed successfully'
    return

if __name__ == '__main__':
    main()
