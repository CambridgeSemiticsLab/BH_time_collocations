from pathlib import Path
repo = Path('/Users/cody/github/CambridgeSemiticsLab/time_collocations')
paths = {
    'time_data': repo.joinpath('results/parsing/time_evaled.json'),
    'phrase_data': repo.joinpath('results/parsing/phrase_parsings.json'),
    'outdir': repo.joinpath('results/analysis'),
}
