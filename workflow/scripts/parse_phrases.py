from tf.fabric import Fabric
from tf.app import use
from parsing.phrase_parser import parse_phrases
from parsing.parse_metrics import examine_parsings

datalocs = {
    'bhsadata': snakemake.input.bhsadata,
    'slot2pos': snakemake.input.slot2pos,
    'plotsdir': snakemake.params.plotsdir,
}

# initialize TF
TF = Fabric(locations=datalocs['bhsadata'], silent='deep')
features = ( 
    'rela code gloss function number '
    'pdp vs vt typ language label st '
    'prs'
)   
API = TF.load(features, silent='deep')
bhsa = use('bhsa', api=API, silent='deep')
bhsa._browse = True
API = bhsa.api

parse_phrases(
    snakemake.input.samples,
    snakemake.output.parsed,
    snakemake.output.notparsed,
    datalocs,
    API=API
)

examine_parsings(
    snakemake.output.parsed,
    snakemake.output.notparsed,
    datalocs,
    snakemake.input.styles,
    snakemake.output.metrics,
    bhsa=bhsa
)
