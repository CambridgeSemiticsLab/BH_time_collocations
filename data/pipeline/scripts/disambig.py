"""
Get disambiguators.
"""

from tf.fabric import Fabric
from disambiguating.accent_sep import get_accent_seps
from disambiguating.para_sep import get_para_seps

paths = {
    "accentseps": snakemake.output.accentseps,
    "paraseps": snakemake.output.paraseps,
}

# initialize TF
TF = Fabric(
    locations=snakemake.input.bhsadata,
    silent='deep'
)
features = ( 
    'pdp ls trailer rela lex mother st'
)   
TF_API = TF.load(features, silent='deep')

# build disambiguations
get_accent_seps(
    paths,
    TF_API
)

get_para_seps(
    paths,
    TF_API
)
