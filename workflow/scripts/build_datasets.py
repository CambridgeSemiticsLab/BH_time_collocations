from dataset.dataset_builder import build_datasets

paths = dict(
    **snakemake.input,
    **snakemake.output,
)

build_datasets(paths)
