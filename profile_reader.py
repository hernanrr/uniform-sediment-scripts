
#!/usr/bin/env python

""" This script reads the profiles stored in the "raw" folder and coverts them
to JSON files in the "input" folder. The JSON file contains a dictionary with
the profile information arranged as a numpy array.  
"""
from __future__ import division
import glob
#import sys
import os
import pdb
import re
import numpy as np
import json
try:
   import cPickle as pickle
except:
   import pickle

class npJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist() # or map(int, obj)
        #if isinstance(obj, datetime.datetime):
        #    return obj.isoformat() # or map(int, obj)
        return json.JSONEncoder.default(self, obj)


# to get the files sorted in the right order
numbers = re.compile(r'(\d+)')
def NumericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


home = os.path.expanduser("~")
sourcepath = (home + '/Documents/Experiments/Data/raw/profiles')
outputpath = (home + '/Documents/Experiments/Data/input/profiles')



def main():
    """Converts the profile info from CSV to JSON files"""
    # Specify working directory as source path and change to it.
    wd = sourcepath
    os.chdir(wd)
    # List all profiles in source path, sorted by feedrate
    profiles = sorted( glob.glob('*.csv'), key=NumericalSort)

    # Create equilibrium and agradation dictionaries:
    eq = {}
    ag = {}
    aggradation = re.compile(r'ag')
    for f in profiles:
        key = f.split('.')[0]
        data_type = ['f8', 'f8', 'f8']
        data_names = 'x, wse, bed'
        # open the file and get the data
        try:
            data = np.genfromtxt(f, dtype = data_type, delimiter=',',
                                 skip_header = 2, names = data_names )
        except ValueError:
            print 'A value error popped up'
            pdb.set_trace()
    # Separate the equilibrium runs from the aggradation runs        
        if aggradation.search(f) is None:
            eq[key] = data
        else: 
            ag[key] = data

    # Dump the dictionaries in their respective places
    # First, the equilibrium runs:
    wd = outputpath + '/equilibrium'
    os.chdir(wd)
    # Create Pickle
    pickle_hdr = 'equilibrium_profiles.pickle'
    with open(pickle_hdr, 'wb') as pickle_outfile:
        pickle.dump(eq, pickle_outfile, -1)    
    
    # Create JSON dump
    # json_header = 'equilibrium_profiles.json'
    # with open(json_header, 'w') as json_outfile:
    #     json.dump(eq, json_outfile, indent = 4, cls=npJSONEncoder)

    # Then the aggradaton runs
    wd = outputpath + '/aggradation'
    os.chdir(wd)

    # Create Pickle
    pickle_hdr = 'aggradation_profiles.pickle'
    with open(pickle_hdr, 'wb') as pickle_outfile:
        pickle.dump(ag, pickle_outfile, -1)    

    # Create JSON dump
    # json_header = 'aggradation_profiles.json'
    # with open(json_header, 'w') as json_outfile:
    #     json.dump(ag, json_outfile, indent = 4, cls=npJSONEncoder)

    print 'Script completed successfully'
    return
    


if __name__ == '__main__':
    main()
