# This script is run on the different nodes when submitting a batch job
# All arguments are set in the Batch_Make_MC_pickles notebook

#imports
import uproot
import glob
import os
import sys
from tqdm import tqdm
import pandas as pd
from pandas import HDFStore, read_hdf
import numpy as np
import collections
from pathlib import Path
import concurrent.futures

# Get parameters defined in Batch_Make_MC_pickles.ipynb
component = sys.argv[1]
isotope = sys.argv[2]
mc_dir = sys.argv[3]
hdf_dir = sys.argv[4]
mc_version = sys.argv[5]
minz = float(sys.argv[6]) #typecast string to float
maxz = float(sys.argv[7]) #typecast string to float
maxr = float(sys.argv[8]) #typecast string to float

_files_dir = mc_dir+component+'_'+isotope+'/'
_hdf_name=hdf_dir+component+'_'+isotope+'.p'


print(component)
print(isotope)
print(mc_dir)
print(hdf_dir)
print(mc_version)

if mc_version!='v2.2.0_G4p10' and mc_version!='v2.4.0_G4p10':
    _endstring='*_g4mc_G4_Sort.root'
else:
    _endstring='*_g4mc_G4p10_Sort.root'

_files = glob.glob(os.path.expanduser(_files_dir) + _endstring)
executor = concurrent.futures.ThreadPoolExecutor(14)

print('Loading metadata.')
meta = {'nevents': [_chunk[1] for _chunk in uproot.iterate(_files, "events/events", ["?"], reportentries=True, entrysteps=100000, 
                                                           outputtype=collections.namedtuple, executor=executor)][-1],
        }

branches= ['Ed','X','Y','Z','type_pri']
print('Extracting branches:', branches)
_chunklist=[]

for _chunk in tqdm(uproot.iterate(_files, "events/events",branches, entrysteps=100000, 
                                       outputtype=pd.DataFrame, executor=executor),total=meta['nevents']/100000):
    _chunk.columns=branches
    _chunk =_chunk[ _chunk['Ed'].apply(lambda x: len(x)==1) ] # Remove events with more than one scatter
    _chunk =_chunk.applymap(lambda x: x[0]) # Convert lists in cells by taking only the 0th element. Use applymap for dataframe
    _chunk.type_pri =_chunk.type_pri.apply(lambda x: x.decode("utf-8")) # Decode bytestrings for primary type. Use apply for dataframe columns
    _chunk['ns']=1
    _chunklist.append(_chunk)
    
print('Concatenating dataframe and saving hdf5 to:', _hdf_name)
df=pd.concat(_chunklist)
hdf = HDFStore(_hdf_name)
hdf.put('data', df)
hdf.close()
print('Done. Exiting')