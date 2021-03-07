"""
This module contains pathways for the primary
data utilized by the project.
"""

from pathlib import Path 

# base paths
home = Path.home()
github = home.joinpath('github')

# project data
project = github.joinpath('CambridgeSemiticsLab/time_collocations')
semvector = project.joinpath('data/semvectors/semvector.pickle')
data = project.joinpath('data')
cxs = data.joinpath('cxs/cxs.pickle')
main_data = data.joinpath('main_dataset')
main_table = main_data.joinpath('dataset.tsv')
annotations = data.joinpath('annotations')
genesis_annotations = annotations.joinpath('Genesis_pilot/genesis_times_annotated.csv')
figs = project.joinpath('results/figures')

# Text-Fabric data
tf_data = {
    'bhsa_app': str(github.joinpath('annotation/app-bhsa')),
    'bhsa': str(github.joinpath('etcbc/bhsa/tf/c')),
    'heads': str(github.joinpath('etcbc/heads/tf/c')),
    'genre': str(github.joinpath('etcbc/genre_synvar/tf/c')),
    'custom': str(project.joinpath('data/bhsa/tf')),
}
