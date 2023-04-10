from pathlib import Path

user_home = Path.home()
REPO_BASE_PATH = 'github/BH_time_collocations/archive/2022-11-02'
REPO_PATH = user_home.joinpath(REPO_BASE_PATH)

paths = {
    'time_dataset': REPO_PATH / 'results/analysis/time_dataset.csv',
    'phrase_dataset': REPO_PATH / 'results/analysis/phrase_dataset.csv',
    'clause_dataset': REPO_PATH / 'results/analysis/clause_dataset.csv',
    'cl_clusters': REPO_PATH / 'results/cl_clusters/clusters.json',
    'time_parses': REPO_PATH / 'results/parsing/time_evaled.json',
    'phrase_parses': REPO_PATH / 'results/parsing/phrase_parsings.json',
    'outdir': REPO_PATH / 'results/analysis',
}
