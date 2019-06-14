# I/O for Geant4 ROOT files
import uproot
import glob
import os
import sys
from tqdm import tqdm_notebook
import pandas as pd
from pandas import HDFStore, read_hdf
import numpy as np
import collections
from pathlib import Path
import concurrent.futures


### Loading ROOT file into data frame
#Using Lutz's UpROOT interface, but scaled down to a single file. It is possible to use multiprocessing if loading several ROOT files from Geant4 with a few modifications.

def g4root_to_df(inputfile='', branches=[], nthreads=2):
    if inputfile == '':
        print('No input file specified. Try again.')
        return 0
    if len(branches) < 1:
        print('No branches specified. Try again.')
        return 0
    else:
        print(inputfile)
        executor = concurrent.futures.ThreadPoolExecutor(2)

        meta = {'nevents': [_chunk[1] for _chunk in uproot.iterate(inputfile, "events/events", ["?"], reportentries=True, entrysteps=1000, 
                                                                   outputtype=collections.namedtuple, executor=executor)][-1],
                }
        
        _chunklist=[]
        for _chunk in tqdm_notebook(uproot.iterate(inputfile, "events/events",branches, entrysteps=1000, 
                                               outputtype=pd.DataFrame, executor=executor),total=meta['nevents']/1000):
            _chunk.columns=branches
            _chunklist.append(_chunk)

        return pd.concat(_chunklist)