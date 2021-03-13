"""
In this module we generate a sample of 
phrases drawn in from the ETCBC's BHSA Hebrew database.
The database is accessed using Text-Fabric.
It is subsequently processed using the scripts and modules indicated below.
"""

from sampling.phrase_sampling import get_phrase_samples

# execute with snakemake arguments
get_phrase_samples(
    snakemake.input.bhsadata,
    snakemake.output.samples,
    snakemake.output.nosamples,
)
