"""
This module generates a custom version of the
BHSA Hebrew database.
"""

from pathlib import Path
from datetime import datetime
from tf.fabric import Fabric
from modify import mod_features
#from time_association import TimeAssociation

home = Path.home()
bhsa_dir = Path(home, 'text-fabric-data/etcbc/bhsa/tf/c')
output_dir = Path('tf_test')

locations = {
    'bhsa': str(bhsa_dir),
    'output': str(output_dir),
}

base_metadata = {
    'modified_by': 'Modifications added by Cody Kingham', 
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
mod_features(locations, base_metadata)
print('*'*20)
msg('done!', 1)

# -- Add Node Features -- 
# -- Find Statistically Heads --
