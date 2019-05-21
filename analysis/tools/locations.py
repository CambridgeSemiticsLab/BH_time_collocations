# Text-Fabric data locations for the analyses

import os.path as path

home = path.expanduser('~')
repo = path.expanduser('~/github/csl/time_collocations')
custom_data = path.join(repo, 'data/tf')
bhsa_data = path.join(home, 'text-fabric-data/etcbc/bhsa/tf/c')
heads_data = path.join(home, 'github/etcbc/heads/tf/c')

data_locations = {'bhsa': bhsa_data,
                  'heads': heads_data,
                  'custom': custom_data}