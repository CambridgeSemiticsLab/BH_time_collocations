# central data locations for the project
import os.path as path
home = path.expanduser('~')
repo = path.expanduser('~/github/CambridgeSemiticsLab/time_collocations')
data = path.join(repo, 'data/')
semvector = path.join(data, 'vectors/semvector.pickle')
cxs = path.join(data, 'cxs/cxs.pickle')
