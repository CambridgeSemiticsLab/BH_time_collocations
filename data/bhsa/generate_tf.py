"""
This module generates a custom version of the
BHSA Hebrew database.
"""

from pathlib import Path
from paths import tf_data
from datetime import datetime
from modify import mod_features
#from time_association import TimeAssociation

base_metadata = {
    'modified': 'Modifications added by Cody Kingham', 
}

start = datetime.now()

def timestamp():
    # give elapsed time 
    return f'{datetime.now()-start}'

def msg(m, indent=0):
    indent = '\t' * indent
    print(f'{indent}{timestamp()}  {m}')

# -- Remap Node Features --
msg('Remapping TF features...')
print('*'*20)
mod_features(tf_data, base_metadata)
print('*'*20)
msg('done!', 1)

# -- Add Node Features -- 
# -- Find Statistically Heads --
