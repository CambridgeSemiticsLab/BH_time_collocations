# Text-Fabric data locations for the project

import os.path as path

home = path.expanduser('~')
repo = path.expanduser('~/github/CambridgeSemiticsLab/time_collocations')
custom_data = path.join(repo, 'data/text_fabric/tf')
bhsa_data = path.join(home, 'text-fabric-data/etcbc/bhsa/tf/c')

data_locations = {'bhsa': bhsa_data,
                  'custom': custom_data}
