"""
Load currently processed CX objects 
from the data directory.
"""

import pickle
from paths import cxs
from .cx import Construction

with open(cxs, 'rb') as infile:
    cxs = pickle.load(infile)
