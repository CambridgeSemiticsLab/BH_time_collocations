from pathlib import Path
repo = Path('/Users/cody/github/CambridgeSemiticsLab/time_collocations')
paths = {
    'time_dataset': repo.joinpath('results/analysis/time_dataset.csv'),
    'phrase_dataset': repo.joinpath('results/analysis/phrase_dataset.csv'),
    'time_parses': repo.joinpath('results/parsing/time_evaled.json'),
    'phrase_parses': repo.joinpath('results/parsing/phrase_parsings.json'),
    'outdir': repo.joinpath('results/analysis'),
}
