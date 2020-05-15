#imports
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

# Define functions to correctly get Run Livetime
import datetime
#import uncertainties.unumpy as unp
import matplotlib as mpl 
from iminuit import Minuit, describe
import datetime
import time
import matplotlib.dates as mdates
import pickle
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.pyplot import cm 

from collections import defaultdict
from collections import OrderedDict

from IPython.display import display, HTML

mpl.rcParams['axes.linewidth'] = 0.75 #set the value globally
mpl.rcParams['xtick.major.width'] = 0.75
mpl.rcParams['ytick.major.width'] = 0.75
mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['pdf.fonttype']=42


mpl.rcParams['figure.figsize'] = (3.49,2.32)
mpl.rcParams['savefig.dpi'] = 300
mpl.rc('font', size= 7)
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)