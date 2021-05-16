from dataset.dataset_builder import build_datasets

paths = {
    'bhsadata': snakemake.input.bhsadata,
    'slot2pos': snakemake.input.slot2pos,
    'timedata': snakemake.input.timedata,
    'tensedata': snakemake.input.tensedata,
    'phrasedata': snakemake.input.phrasedata,
    'functions': snakemake.input.functions,
    'timedataset': snakemake.output.timedataset,
    'phrasedataset': snakemake.output.phrasedataset,
}

build_datasets(paths)
