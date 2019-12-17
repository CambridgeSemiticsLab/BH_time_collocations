import collections
from tf.fabric import Fabric
from tf.compose import modify
from weqetal import convert_tense

def mod_features(locations, metadata):
    """Remap node features in BHSA

    BHSA contains features of  nodes that 
    are sometimes inaccurate or not quite
    fit for purpose for this study. In this
    module, we modify the dataset with Text-Fabric
    and export the new version to a new directory.
    """

    # setup data and methods
    bhsa, output = locations['bhsa'], locations['custom']
    TF = Fabric(locations=bhsa, silent=True)
    api = TF.load('vt mother pdp lex')

    # list node features to be modified
    to_modify = {
        'nodeFeatures': {
            'function': {
                849296:'Loca', 
                825329:'Loca',
                828081:'Cmpl',
                774349:'Adju',
                774352:'Adju',
                775948:'Adju',
                775985:'Adju',
                876172:'Adju',
                881665:'Objc', # phrase belongs with previous as adjectival element
                870273:'Time', # prep and conj belong with time phrase
                870274:'Time', # modifier "KBR" belongs with time phrase
            },
           'vt': {
            }
        },
    }

    # modify verb tenses to add weqetal
    for verb in api.F.pdp.s('verb'):
        to_modify['nodeFeatures']['vt'][verb] = convert_tense(verb, api)

    # enact the changes
    modify(
        bhsa,
        output,
        addFeatures=to_modify,
        featureMeta={feat:metadata for feat in ('function','vt')},
        deleteFeatures=True, # export only modified features
    )
