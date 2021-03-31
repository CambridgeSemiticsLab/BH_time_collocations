from tf.fabric import Fabric
from tf.app import use
from parsing.time_parser import parse_times
from parsing.time_metrics import examine_times

datalocs = {
    'bhsadata': snakemake.input.bhsadata,
    'plotsdir': snakemake.params.plotsdir,
    'ph_parsings': snakemake.input.ph_parses,
}

# initialize TF
TF = Fabric(locations=datalocs['bhsadata'], silent='deep')
features = ( 
    'rela code gloss function number '
    'pdp vs vt typ language label st '
    'prs nu'
)   
API = TF.load(features, silent='deep')
bhsa = use('bhsa', api=API, silent='deep')
bhsa._browse = True # ensures API.pretty outputs HTML strings
API = bhsa.api

parse_times(
    snakemake.input.ph_parses,
    snakemake.output.parsed,
    snakemake.output.notparsed,
    datalocs,
    API=API
)

examine_times(
    snakemake.output.parsed,
    snakemake.output.notparsed,
    datalocs,
    snakemake.input.styles,
    snakemake.output.metrics,
    bhsa=bhsa
)
