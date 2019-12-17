"""
This module contains pathways for the primary
data utilized by the project.
"""

from pathlib import Path 

# base paths
home = path.home()
github = home.joinpath('gitub')

# project data
project = github.joinpath('CambridgeSemiticsLab/time_collocations')
semvector = project.joinpath('data/vectors/semvector.pickle')
cxs = project.joinpath('data/cxs/cxs.pickle')

# Text-Fabric Data
tf_data = {
    'bhsa': github.joinpath('etcbc/bhsa/tf/c')
    'heads': github.joinpath('etcbc/heads/tf/c')
    'custom': project.joinpath('data/bhsa/tf')
}
