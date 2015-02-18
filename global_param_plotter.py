#!/usr/bin/env python
""" This script plots the global parameters

"""
from __future__ import division
import glob
import sys
import os
import pdb
import re 
import numpy as np
try:
    import cPickle as pickle
except:
    import pickle
import itertools

# Third party modules
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Rectangle
from matplotlib import rc
from matplotlib import rcParams
rc('font', size = 10, **{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rcParams['pdf.fonttype'] = 42
rcParams['text.color'] = 'black'
rcParams['text.latex.preamble'] = [
    r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

numbers = re.compile(r'(\d+)')
def NumericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
# Some Constants

# g = 9.81 # Gravity,[ m / s ** 2 ]
# nu =  1.0e-6 # Dynamic viscosity of water [ m **2 / s ]
# rho = 1000 # Density of water [ kg / m ** 3 ]

# # Ashida Michiue bedload transport relation.
# alpha = 17
# tau_star_c = 0.05

# # Meyer-Peter and Muller, corrected by Wong bedload transport relation
# alpha_MPM = 3.97
# tau_star_c_MPM = 0.05
# n_MPM = 1.5

# #Flume stuff
# B0 = 0.19 # Meters of channel width
# rho_s = 2650 # kg/m3
# R = 1.65
# D = 1.11/1000 # meters, characteristic grain size



# # Plot the data
# fig = plt.figure(figsize=(6,6), tight_layout=True)
# ax1 = fig.add_subplot(211, aspect = 'equal')
# for key in d.keys():
#     data = d[key]['tau_star_bs'] #- tau_star_c
#     x = excess_shields[key]
#     y = q_b_star[key]
#     xmin = np.amin(data)
#     xmax = np.amax(data)
#     height = 0.001
#     print key, xmin, xmax
#     ax1.add_patch(Rectangle((xmin,y-height/2), (xmax-xmin), height, facecolor
#                   = 'blue', edgecolor = 'blue', alpha=0.5 ))
#     ax1.scatter(x, y, s=48, facecolor = 'blue', edgecolors='grey', lw = 0)
# #    ax1.axhline(y = y, xmin=xmin, xmax=xmax, linewidth=4, color='blue',
# #                alpha=0.5)
#     ax1.set_xscale('log')
#     ax1.set_yscale('log')



# ax2 = fig.add_subplot(212, aspect = 'equal')
# for key in d.keys():
#     data = d[key]['tau_star_bs'] #- tau_star_c
#     x = total_shields[key]
#     y = phi[key]
#     ax2.scatter(x, y, s=48, facecolor = 'blue', edgecolors='grey', lw = 0)








# #title1 = 'Sonar measurements' 
# #fig.suptitle(title1, fontsize=24)





# Establish bedload relations


def load_pickle(sourcepath, run, suffix):
    """Load the profiles from the pickle"""
    wd = os.path.join(sourcepath, run)
    os.chdir(wd)
    fpickle = '{}_{}.pickle'.format(run, suffix)
    with open(fpickle, 'rb') as infile:
        pkl = pickle.load(infile)
    return pkl 

def save_fig(fig, outputpath, run, figname):
    """Saves the figure to disk"""
    wd = os.path.join(outputpath, run, 'plots', 'global_parameters')
    os.chdir(wd)
    f = '{}_{}.pdf'.format(run, figname)
    fig.savefig(f, dpi=300, format='pdf', transparent=True,
                bbox_inches='tight', pad_inches=0.1, frameon=False)
    print 'File : {} written to {}'.format(f, wd)

def ashida_michiue(tau, tau_star_c=0.05):
    """Computes the Ashida & Michiue Bedload"""
    qbstar = 17 * ( tau - tau_star_c) * (np.sqrt(tau) - np.sqrt(tau_star_c) )
    return qbstar

def mpm_wong(tau, tau_star_c=0.0495):
    """Computes the Ashida & Michiue Bedload"""
    qbstar = 3.97 * ( tau - tau_star_c) ** 1.5
    return qbstar


am_eq = (r'$q^{*}_{bs} = \num{17} (\tau^{*}_{\text{bs}} - \num{0.05})'+\
         r'(\sqrt{\tau^{*}_{\text{bs}}} - \sqrt{\num{0.05}} )$')
mpm_eq = r'$q^{*}_{bs} = \num{3.97} (\tau^{*}_{\text{bs}} - \num{0.0495})^{\num{1.5}}$'

bedload_rels = {}
bedload_rels['A-M'] = {'alpha':np.float(17),
                       'tau_star_c':np.float(0.05),
                       'eq': am_eq,
                       'name': 'Ashida & Michiue',
                       'fcall': ashida_michiue,
                       'color': 'black',
                       'linestyle': '--',
                   }
bedload_rels['MPM-W'] = {'alpha':np.float(3.97),
                         'tau_star_c':np.float(0.0495),
                         'n':np.float(1.5),
                         'eq': mpm_eq, 
                         'name': 'MPM, modified by Wong',
                         'fcall':mpm_wong,
                         'color': 'black',
                         'linestyle': '-',
                     }

def plot_load_relations(ax, bedload_rels):
    """Plots load relation curves"""
    # Specify plotting domain
    tau = np.linspace(0.05, .6, 500)
    for key in bedload_rels:
        qbstar = bedload_rels[key]['fcall'](tau)
        label = r'{}'.format( bedload_rels[key]['eq'] )
        color = bedload_rels[key]['color']
        linestyle = bedload_rels[key]['linestyle']
        ax.loglog(tau, qbstar, label=label, color = color, linestyle =
                  linestyle)
    return

def plot_stresses(ax, gp, stats):
    """Plot the experimental points' stresses"""
    for key in sorted(stats, key=NumericalSort):
        x = np.mean(stats[key]['Mean']['taub_star_s'])
        y = stats[key]['Meta']['qb_star']
        t = 0.001 # thickness
        units = r'\,\si{\g \per \minute}'
        label = r'\num{} {}'.format(key, units)
        xmin = np.amin(stats[key]['Mean']['taub_star_s'])
        xmax = np.amax(stats[key]['Mean']['taub_star_s'])
        ax.add_patch(Rectangle((xmin, y-t/2), (xmax-xmin), t, alpha=0.5))
        ax.plot(x, y, markersize=6, label=label, linestyle = 'None',
                marker=u'o', markeredgecolor='gray', markeredgewidth=0.2)
#        ax.scatter(x, y, s=48, facecolor = 'blue', edgecolor = 'blue', lw = 0)
    return


def label_axes(ax, xlabel, ylabel, fontsize):
    """Label axes and title figure"""
    ax.set_xlabel(xlabel, fontsize = fontsize)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    return

def format_plot(fig, xlim, ylim, xscale, yscale):
    """Prettify plots"""
    spines_to_remove = ['top', 'right']
    spines_to_keep = ['bottom', 'left']
    for ax in fig.axes:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.xaxis.set_visible(True)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        for spine in spines_to_remove:
            ax.spines[spine].set_visible(False)
        for spine in spines_to_keep:
            ax.spines[spine].set_color('black')
        ax.tick_params(\
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            left='on',         # ticks along the left edge are off
            right='off',         # ticks along the right edge are off
            labelbottom='on', # labels along the bottom edge are off
            labelleft='on',  # labels along the bottom edge are off
            labelcolor = 'black') # Labels are black
    #for ax in [ax.xaxis, ax.yaxis]:
    #    ax.set_major_formatter(ScalarFormatter())

    return
    

def plot_phi(ax, stats):
    """Plot the experimental points' stresses ratio"""
    for key in sorted(stats, key=NumericalSort):
        x = np.mean(stats[key]['Mean']['taub_star'])
        y = np.mean(stats[key]['Mean']['phi'])
        height = np.ptp(stats[key]['Mean']['phi'])
        w = 0.001 # width
        ymin = np.amin(stats[key]['Mean']['phi'])
        ymax = np.amax(stats[key]['Mean']['phi'])
        units = r'\,\si{\g \per \minute}'
        label = r'\num{} {}'.format(key, units)
        ax.add_patch(Rectangle((x-w/2, ymin), w, height, alpha=0.5))
        ax.plot(x, y, markersize=6, label=label, linestyle = 'None',
                marker=u'o', markeredgecolor='gray', markeredgewidth=0.2)

    return


def plot_phi_average(ax, stats):
    """Plot the average experimental points' stresses ratio"""
    for key in sorted(stats, key=NumericalSort):
        x = np.mean(stats[key]['Mean']['taub_star'])
        y = ( np.mean(stats[key]['Mean']['taub_star_s']) 
              / np.mean(stats[key]['Mean']['taub_star']) )
        ax.scatter(x, y, s=48, lw = 0, label=key)
    return

def plot_friction(ax):
    """Plot the friction relation"""
    alpha_r, nk = np.float(8.1), np.float(1/6)
    rb_slash_ks = np.append(np.linspace(1, 10, 9, endpoint = False), 
                            np.linspace(10, 100, 10, endpoint = True))
    Cz = alpha_r * rb_slash_ks ** nk
    Cfe = 1 / ( Cz ** 2 )
    label = (r'$C_{\text{f}_{\text{bs}}} = \alpha_r \left(' + \
             r'\frac{R_{\text{H}_{\text{bs}}}}{k_{\text{s}}} \right)^{1/6}$')
    ax.plot(rb_slash_ks, Cfe, label = label)
    return
    
def load_points(pointspath):
    """Load the experimental points of MPM, Gilbert and Viparelli"""
    # Change directory to location of points
    os.chdir(pointspath)
    f_mpm = 'MPM_points.csv'
    f_gilbert = 'gilbert_points.csv'
    f_viparelli = 'viparelli_points.csv'

    mpm = np.loadtxt(f_mpm, skiprows = 1, delimiter=',')
    gilbert = np.loadtxt(f_gilbert, skiprows = 1, delimiter=',')
    viparelli = np.loadtxt(f_viparelli, skiprows = 1, delimiter=',')
    return mpm, gilbert, viparelli

def recompute_rb_over_ks(Rb_over_ks, nk=1.5):
    """Recomputes Rb_over_ks. Not sure what's going on here.  This seems necessary
    for the MPM and Gilbert sets. Nothing to worry about

    """
    return Rb_over_ks * 2 / nk
    

def compute_chezy_friction_coeff(Cz):
    """Computes the Chezy friction coefficient given the Chezy resistance 
    factor.

    """
    return 1./ ( Cz ** 2 )  

def process_points(mpm, gilbert, viparelli):
    """Process the points to make them suitable for plotting"""
    for points in [mpm, gilbert]:
        points[:,0] = recompute_rb_over_ks(points[:,0])
    for points in [mpm, gilbert, viparelli]:
        points[:,1] = compute_chezy_friction_coeff(points[:,1])
    return mpm, gilbert, viparelli
    
    

def plot_others_points(ax, x, y, label, marker='o'):
    """Process and plot experimental points"""
    ax.plot(x, y, label=label, linestyle = 'None', marker = marker,
    markersize=6, markeredgecolor = 'gray', markeredgewidth=0.2)
    return
    
def plot_my_points(ax, s):
    """Plots my points"""
    for key in sorted(s, key=NumericalSort):
        label = r'Our data'
        x = np.mean(s[key]['Mean']['Rhb_s']) / s[key]['Meta']['ks']
        y = np.mean(s[key]['Mean']['Cfbs'])
        ax.plot(x, y, label=label, linestyle = 'None', marker = u'o',
                markersize=6, markerfacecolor='#8da0cb',
                markeredgecolor='gray', markeredgewidth=0.2)
    return

def print_table(gp, stats, run):
    """Write line to file"""
    # Define LaTeX headers
    c1 = r'G_f/(\si{\g \per \min}),'
    c2 = r'q_b^*,'
    c3 = r'u^*_{\text{b}}/(\si{\m\per\s}),'
    c4 = r'u^*_{\text{bs}}/(\si{\m\per\s}),'
    c5= r'$\tau^*_{\text{b}}$,'
    c6 = r'$\tau^*_{\text{bs}}$,'
    c7 = r'$\phi$,'
    c8 = r'H / \si{\cm},'
    c9 = r'U / \si{\m \per \s},'
    c10 = r'Fr'
    hdr = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10
    table=[]
    for key in sorted(stats, key=NumericalSort):
        Gf = np.int(key)
        qbstar = stats[key]['Meta']['qb_star']
        ub_star = np.mean(stats[key]['Mean']['ub_star']) *100
        ub_star_s = np.mean(stats[key]['Mean']['ub_star_s']) * 100
        tau_star_b = np.mean(stats[key]['Mean']['taub_star']) 
        tau_star_bs = np.mean(stats[key]['Mean']['taub_star_s'])
        phi = np.mean(stats[key]['Mean']['phi'])
        H = np.mean(stats[key]['Mean']['H']) * 100.
        U = np.mean(stats[key]['Mean']['U'])
        Fr = np.mean(stats[key]['Mean']['Fr'])
        row = np.array([Gf, qbstar, ub_star, ub_star_s, tau_star_b,
                        tau_star_bs, phi, H, U, Fr])
        table.append(row)
    output = np.vstack(table)
    fname = '/Users/ricardo/Documents/Experiments/Data/output/profiles/equilibrium/tables/table_eq.csv'
    np.savetxt(fname, output, delimiter=',', header=hdr)
    return
    

def main():
    """Main routine for plotting global parameters"""
    print 'Script started' 
    home = os.path.expanduser("~")
    sourcepath = (home + '/Documents/Experiments/Data/processed/profiles')
    outputpath = (home + '/Documents/Experiments/Data/output/profiles')
    pointspath = (home + '/Documents/Experiments/Data/scripts')
    runs = ['equilibrium', 'aggradation']
    for run in runs:
        # Load the pickle with the relevant global parameters
        gp = load_pickle(sourcepath, run, 'global_parameters')
        # Load the pickle with the stats summary
        stats = load_pickle(sourcepath, run, 'global_stats_summary')
        # Print equilibrium results table to file
        print_table(gp, stats, run)
        print 'Plotting experimental points'
        # Create a figure for the plots of experimental points
        fig = plt.figure(tight_layout=True)
        # Create axis for existing load relation curves
        ax1 = fig.add_subplot(111, aspect = 'equal')
        # Plot the load relation curves (Ashida & Michiue and MPM-W, for now)
        plot_load_relations(ax1, bedload_rels)
        # Plot the relevant points
        plot_stresses(ax1, gp, stats)
        # Add title and labels
        xlabel1 = r'$\tau_{bs}^*$'# - \tau_{ref}^*$'
        ylabel1 = r'$q{_b}{^*}$'#r'q$^*_b$'
        label_axes(ax1, xlabel1, ylabel1, 10)
        #title1 = r''
        #fig.suptitle(title1, fontsize = 20)
        ax1.legend(fontsize=8, loc='upper left', numpoints = 1, frameon=False)
        format_plot(fig, (0, 1), (0.05, 10), 'log', 'log')

        # Save the figure of experimental points
        save_fig(fig, outputpath, run, 'experimental_points')

###############
# Plot of phi #
###############
        print 'Plotting phi'
        # Create a plot for phi
        phi = plt.figure()
        ax2 = phi.add_subplot(111)
        # Plot phi
        plot_phi(ax2, stats)
        ax2.legend(fontsize=8, loc='upper left', numpoints = 1, frameon=False)
        # Label the phi plot
        xlabel2 = r'$\tau_{b}^*$'# - \tau_{ref}^*$'
        ylabel2 = r'$\varphi$'#r'q$^*_b$'
        label_axes(ax2, xlabel2, ylabel2, 10)
        format_plot(phi, (0, 1), (0, 1), 'linear', 'linear')
        
        # Save the figure of phi
        save_fig(phi, outputpath, run, 'skin_friction_fraction')


# ##########################
# # plot of phi of average #
# ##########################
#         print 'Plotting phi'
#         # Create a plot for phi
#         phiavg = plt.figure()
#         ax3 = phiavg.add_subplot(111)
#         # Plot phi
#         plot_phi_average(ax3, stats)
#         # Label the phi plot
#         xlabel3 = r'$\tau_{b}^*$'# - \tau_{ref}^*$'
#         ylabel3 = r'$\varphi$'#r'q$^*_b$'
#         label_axes(ax3, xlabel3, ylabel3, 10)
#         format_plot(phiavg, (0, 1), (0, 1), 'linear', 'linear')
        
#         # Save the figure of phi
#         save_fig(phiavg, outputpath, run, 'skin_friction_fraction_avg')
#         #Close all figures


#################
# plot friction #
#################
        print 'Plotting friction'
        # Create a plot for the friction
        friction = plt.figure()
        ax4 = friction.add_subplot(111)
        # Plot friction relation
        plot_friction(ax4)
        # Get the experimental points from MPM, Gilbert and Viparelli
        points = load_points(pointspath)
        # Process others' points
        mpm, gilbert, viparelli = process_points(*points)
        # Plot MPM
        plot_others_points(ax4, mpm[:,0], mpm[:,1], 'MPM', marker='v')
        # Plot Gilbert
        plot_others_points(ax4, gilbert[:,0], gilbert[:,1], 'Gilbert',
                           marker='x')
        # Plot Viparelli
        plot_others_points(ax4, viparelli[:,0], viparelli[:,1], 'Viparelli',
                           marker='s')
        # Plot my points
        #
        #
        plot_my_points(ax4, stats)
        # http://stackoverflow.com/questions/13303928/how-to-make-custom
        #-legend-in-matplotlib
        #Get artists and labels for legend and chose which ones to display
        handles, labels = ax4.get_legend_handles_labels()
        display = (0, 1, 2, 3, 4)
        #Create legend from custom artist/label lists
        ax4.legend([handle for i ,handle in enumerate(handles) if i in display],
                  [label for i, label in enumerate(labels) if i in display],
                   loc='upper left', frameon=False, fontsize=8 )

        
        #ax4.legend(fontsize=8, loc='upper left', numpoints = 1, frameon=False)
        xlabel4 = r'$R_{\text{H}_{\text{bs}}} / k_{\text{s}}$'
        ylabel4 = r'$C_{\text{f}_{\text{bs}}}$'
        label_axes(ax4, xlabel4, ylabel4, 10)
        format_plot(friction, (1, 100), (1e-3, 1e-1), 'log', 'log')
        # Save the figure of friction
        save_fig(friction, outputpath, run, 'friction_plots')

        plt.close('all')
    
    print 'Script completed successfully'
    return


if __name__ == '__main__':
    main()
 
