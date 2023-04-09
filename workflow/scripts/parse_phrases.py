from tf.fabric import Fabric
from tf.app import use
from parsing.phatom_parser import parse_phrases
from parsing.phatom_composer import compose_phrases
from parsing.parse_metrics import examine_parsings

paths = {
    'samples': snakemake.input.samples,
    'bhsadata': snakemake.input.bhsadata,
    'slot2pos': snakemake.input.slot2pos,
    'paraseps': snakemake.input.paraseps,
    'editedges': snakemake.input.editedges,
    'styles': snakemake.input.styles,
    'plotsdir': snakemake.params.plotsdir,
    'parsed_atoms': snakemake.output.parsed_atoms,
    'unparsed_atoms': snakemake.output.unparsed_atoms,
    'parsed_phrases': snakemake.output.parsed_phrases,
    'parsedmets': snakemake.output.parsemetrics,
    'notparsedmets': snakemake.output.notparsemetrics,
}

# initialize TF
TF = Fabric(locations=paths['bhsadata'], silent='deep')
features = ( 
    'rela code gloss function number '
    'pdp vs vt typ language label st '
    'prs ls mother'
)   
API = TF.load(features, silent='deep')
bhsa = use('bhsa', api=API, silent='deep')
bhsa._browse = True
API = bhsa.api

# run the parser, parse composer, and parse metrics
print('Parsing phrase_atoms...')
parse_phrases(paths, API)

print('Composing phrase_atoms...')
compose_phrases(paths, API)

print('Examining final parses...')
examine_parsings(paths, bhsa)
