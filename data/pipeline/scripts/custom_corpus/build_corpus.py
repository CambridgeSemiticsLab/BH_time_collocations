"""Module for constructing a custom corpus from ETCBC's BHSA."""

import yaml

from pathlib import Path

from kingham_thesis.data_pipeline.corpus_building.corpus_builder import ThesisCorpusBuilder

from config import THESIS_CORPUS_PARAMS


def load_yaml(filepath):
    """Load yaml config as dict."""
    with open(filepath, 'r') as infile:
        return yaml.load(infile, Loader=yaml.FullLoader)


# load BHSA generic metadata
BHSA_YAML = snakemake.params.bhsa_repo / 'yaml'
BHSA_GENERIC = BHSA_YAML / 'generic.yaml'
BHSA_GENERIC_META = load_yaml(BHSA_GENERIC)
BHSA_GENERIC_META.update({
    'dateWritten': None,
    'writtenBy': None,
})


DATA_LOCATIONS = [
    snakemake.params.bhsa_data_path,
    snakemake.params.bhsa_genre_path,
    snakemake.params.bsha_heads_path,
]

annotation_files = [Path(file) for file in snakemake.input.annotations]

corpus_builder = ThesisCorpusBuilder(
    locations=DATA_LOCATIONS,
    filter_metadata=BHSA_GENERIC_META,
    annotation_files=annotation_files,
    **THESIS_CORPUS_PARAMS,
)

corpus_builder.build(str(snakemake.output))
