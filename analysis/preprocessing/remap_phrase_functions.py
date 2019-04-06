'''
-- Correcting Time Functions in BHSA -- 
In this script, phrase functions are redrawn
based on the manual corrections stored in the 
dictionary named `newfunctions`
'''

import collections
from tf.fabric import Fabric

dirs = {'bhsa': '~/text-fabric-data/etcbc/bhsa/tf/c', 
        'export': '~/github/csl/time_collocations/data'}

# load TF api + functions
TF = Fabric(locations=dirs['bhsa'], silent=True)
api = TF.load('''

function

''', silent=True)

# draw phrase functions from BHSA base data
phrase2funct = dict((ph, api.F.function.v(ph)) for ph in api.F.otype.s('phrase'))
note = {}

# map corrections/changes here
newfunctions = {849296:'Loca',
                825329:'Loca',
                828081:'Cmpl',
                774349:'Adju',
                774352:'Adju',
                775948:'Adju',
                775985:'Adju',
                876172:'Adju',}

phrase2funct.update(newfunctions)
for phrase, new_funct in newfunctions.items():
    note[phrase] = f'function changed from Time to {new_funct}'

# export
meta = {'':{'source': 'https://github.com/etcbc/bhsa',
            'origin': 'Made by the ETCBC of the Vrije Universiteit Amsterdam; edited by Cody Kingham'},
        'note':{'valueType':'str',
                'description':'notes on objects for tracking issues throughout my research'},
        'function':{'valueType':'str'}}

TFexport = Fabric(locations=dirs['export'], silent=True)
TFexport.save(metaData=meta, nodeFeatures={'function':phrase2funct, 'note':note})

print('remapped functions ready...')