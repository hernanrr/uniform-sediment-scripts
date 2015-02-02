#!/usr/bin/env python

""" This script reads the profiles stored in the "raw" folder and coverts them
to JSON files in the "input" folder. The JSON file contains a dictionary with
the profile information arranged as a numpy array.  
"""
from __future__ import division
import glob
import os
import pdb
import re
import numpy as np
import json
try:
   import cPickle as pickle
except:
   import pickle

# Set input and output folder locations
home = os.path.expanduser("~")
sourcepath = (home + '/Documents/Experiments/Data/raw/flume')
outputpath = (home + '/Documents/Experiments/Data/input/flume')



def main():
    """Converts the profile info from CSV to pickles"""
    # Specify working directory as source path and change to it.
    wd = sourcepath
    os.chdir(wd)
    # List all profiles in source path, sorted by feedrate
    profile = 'flume.csv'
    pdb.set_trace()
    data_type = ['f8', 'f8', 'f8']
    data_names = 'x, wse, bed'
    # open the file and get the data
    try:
       data = np.genfromtxt(profile, dtype = data_type, delimiter=',',
                            skip_header = 1, names = data_names )
    except ValueError:
       print 'A value error popped up'
       pdb.set_trace()

    # Dump the profile in the corresponding input folder
    wd = outputpath
    os.chdir(wd)

    # Create Pickle
    pickle_hdr = 'flume_profile.pickle'
    with open(pickle_hdr, 'wb') as pickle_outfile:
        pickle.dump(data, pickle_outfile, -1)    

    print 'Script completed successfully'
    pdb.set_trace()
    return
    


if __name__ == '__main__':
    main()
