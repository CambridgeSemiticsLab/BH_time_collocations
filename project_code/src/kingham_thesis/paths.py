"""Manages pathing for the project."""

from pathlib import Path

user_home = Path.home()
REPO_NAME = 'BH_time_collocations'
repo = user_home.joinpath('github') / REPO_NAME
data = repo / 'data'
parsing = data / 'parsing'
tables = data / 'tables'

paths = {
    'time_dataset': (tables / 'time_dataset.csv'),
    'phrase_dataset': (tables / 'phrase_dataset.csv'),
    'clause_dataset': (tables / 'clause_dataset.csv'),
    'cl_clusters': (data / 'cl_clusters/clusters.json'),
    'time_parses': (parsing / 'time_evaled.json'),
    'phrase_parses': (parsing / 'phrase_parsings.json'),
}

