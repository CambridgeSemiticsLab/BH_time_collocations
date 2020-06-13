'''
Sets up an instance of TF along by 
loading standard features + custom ones.
'''

from paths import tf_data
from tf.app import use
from tf.fabric import Fabric

standard = ''' 
pdp vs vt
lex language gloss voc_lex voc_lex_utf8
function typ rela code
number label book
'''
plus = '''
function2
lex_sbl lex_sbl_l g_cons_utf8 trailer_utf8
lex_utf8
nu sp prs ls st
mother
freq_lex
kind 
'''
default = standard + plus

locations = list(tf_data.values())[1:]

def load_tf(features='', silent=False, **kwargs):
    """Load an instance of TF with standard and other features
    
    Args:
        features: str of features to load
        silent: str or bool for tf.Fabric load report
        **kwargs: kwargs to pass to tf.app.use

    Returns:
        3-tuple of (TF, api, A)
    """
    TF = Fabric(locations=locations, silent=silent)
    features = default + features
    api = TF.load(features, silent=silent)
    A = use(tf_data['bhsa_app'], api=api, silent=True, **kwargs)
    return (TF, api, A)
