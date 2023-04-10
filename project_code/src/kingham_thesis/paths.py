from pathlib import Path

user_home = Path.home()
REPO_NAME = 'BH_time_collocations'
repo = user_home.joinpath('github') / REPO_NAME

paths = {
    'time_dataset': repo.joinpath('results/analysis/time_dataset.csv'),
    'phrase_dataset': repo.joinpath('results/analysis/phrase_dataset.csv'),
    'clause_dataset': repo.joinpath('results/analysis/clause_dataset.csv'),
    'cl_clusters': repo.joinpath('results/cl_clusters/clusters.json'),
    'time_parses': repo.joinpath('results/parsing/time_evaled.json'),
    'phrase_parses': repo.joinpath('results/parsing/phrase_parsings.json'),
    'outdir': repo.joinpath('results/analysis'),
}
