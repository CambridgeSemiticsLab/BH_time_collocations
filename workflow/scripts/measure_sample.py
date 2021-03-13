from sampling.metrics import measure_phrases

measure_phrases(
    snakemake.input.samples,
    snakemake.input.nosamples,
    snakemake.input.bhsadata,
    snakemake.output.metrics,
    snakemake.params.outdir
)
