"""
This module generates a custom version of the
BHSA Hebrew database.
"""

import sys
import os
import glob
from paths import tf_data
from datetime import datetime
from modify import mod_features
from add import to_graph
#from time_association import TimeAssociation

base_metadata = {
    'source': 'https://github.com/etcbc/bhsa',
    'data_version': 'BHSA version c',
    'origin': 'Made by the ETCBC of the Vrije Universiteit Amsterdam',
    'modified': 'Modified by Cody Kingham for the PhD project, Time Collocations',
}

start = datetime.now()

def timestamp():
    # give elapsed time 
    return f'{datetime.now()-start}'

def msg(m, indent=0):
    indent = ' | \t' * indent
    print(f'{indent}{timestamp()} {m}')

do_assoc = '-assoc' in sys.argv
assoc_files = {os.path.join(tf_data['custom'], f) for f in ('funct_assoc.tf','top_assoc.tf')}

# clean out old data
out_dir = os.path.join(tf_data['custom'], '*.tf')
for file in glob.glob(out_dir):
    if file in assoc_files and not do_assoc:
        continue
    os.remove(file)

# -- Remap Node Features --
print()
msg('Editing node features...')
mod_features(tf_data, base_metadata, do_assoc=do_assoc)
msg('done', 1)

msg('Adding new nodes and edges...')
to_graph(tf_data, base_metadata)
msg('done', 1)

msg('Modifications complete!')
print()
