from dataset.dataset_builder import build_dataset

paths = {
    'bhsadata': snakemake.input.bhsadata,
    'timedata': snakemake.input.timedata,
    'phrasedata': snakemake.input.phrasedata,
    'functions': snakemake.input.functions,
    'timedataset': snakemake.output.timedataset,
    'phrasedataset': snakemake.output.phrasedataset,
}

build_dataset(paths)
