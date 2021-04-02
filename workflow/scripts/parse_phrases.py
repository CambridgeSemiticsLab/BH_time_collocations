from tf.fabric import Fabric
from tf.app import use
from parsing.phrase_parser import parse_phrases
from parsing.parse_metrics import examine_parsings

paths = {
    'samples': snakemake.input.samples,
    'bhsadata': snakemake.input.bhsadata,
    'slot2pos': snakemake.input.slot2pos,
    'paraseps': snakemake.input.paraseps,
    'styles': snakemake.input.styles,
    'plotsdir': snakemake.params.plotsdir,
    'parsed': snakemake.output.parsed,
    'notparsed': snakemake.output.notparsed,
    'parsedmets': snakemake.output.parsemetrics,
    'notparsedmets': snakemake.output.notparsemetrics,
}

# initialize TF
TF = Fabric(locations=paths['bhsadata'], silent='deep')
features = ( 
    'rela code gloss function number '
    'pdp vs vt typ language label st '
    'prs ls'
)   
API = TF.load(features, silent='deep')
bhsa = use('bhsa', api=API, silent='deep')
bhsa._browse = True
API = bhsa.api

# run the parser and parse metrics
parse_phrases(paths, API)
examine_parsings(paths, bhsa)
