from tf.fabric import Fabric
from tf.app import use
from parsing.time_parser import parse_times
from parsing.time_metrics import examine_times

paths = { 
    'bhsadata': snakemake.input.bhsadata,
    'ph_parses': snakemake.input.ph_parses,
    'styles': snakemake.input.styles,
    'plotsdir': snakemake.params.plotsdir,
    'parsed': snakemake.output.parsed,    
    'notparsed': snakemake.output.notparsed,
    'metrics': snakemake.output.parsemets,
    'catmetrics': snakemake.output.catmets,
}

# initialize TF
TF = Fabric(locations=paths['bhsadata'], silent='deep')
features = ( 
    'rela code gloss function number '
    'pdp vs vt typ language label st '
    'prs nu mother'
)   
API = TF.load(features, silent='deep')
bhsa = use('bhsa', api=API, silent='deep')
bhsa._browse = True # ensures API.pretty outputs HTML strings
API = bhsa.api

parse_times(
    paths,
    API=API
)

examine_times(
    paths,
    bhsa=bhsa
)
