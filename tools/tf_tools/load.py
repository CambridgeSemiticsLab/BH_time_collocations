'''
Sets up an instance of TF along by 
loading standard features + custom ones.
'''

import os.path as path
from tf.app import use
from tf.fabric import Fabric

home = path.expanduser('~')
repo = path.expanduser('~/github/CambridgeSemiticsLab/time_collocations')
custom_data = path.join(repo, 'data/text_fabric/tf')
bhsa_data = path.join(home, 'text-fabric-data/etcbc/bhsa/tf/c')
tf_data = [bhsa_data, custom_data]

standard = ''' 
pdp vs vt
lex language gloss voc_lex voc_lex_utf8
function typ rela
number label book
'''
plus = '''
g_cons_utf8 trailer_utf8
nu sp prs ls st
mother
freq_lex
kind 
top_assoc
'''
default = standard + plus

def load_tf(features='', silent=False):
    """Load an instance of TF with standard and other features"""
    TF = Fabric(locations=tf_data, silent=silent)
    features = default + features
    api = TF.load(features, silent=silent)
    A = use('bhsa', api=api, silent=True)
    return (TF, api, A)
