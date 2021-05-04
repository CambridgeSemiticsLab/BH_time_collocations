from dataset.dataset_builder import build_dataset

paths = {
    'bhsadata': snakemake.input.bhsadata,
    'timedata': snakemake.input.timedata,
    'phrasedata': snakemake.input.phrasedata,
    'dataset': snakemake.output.dataset,
}

build_dataset(paths)
