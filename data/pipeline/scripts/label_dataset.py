"""This module labels the analysis dataset."""

from tf.fabric import Fabric

from labeling.autolabeler import AutoLabeler
from labeling.query_labeler import QueryLabeler
from labeling.params import (
    annotation_obj_specs, label_specs, label_queries,
)
from labeling.annotation_sheets import BasicAnnotationSheet

# configure resources
tf_fabric = Fabric(
    locations=snakemake.input.corpus,
    silent="deep",
)
tf_api = tf_fabric.loadAll()

processors = [
    QueryLabeler(tf_fabric, label_queries),
]


# get labels
labeler = AutoLabeler(
    outdir=snakemake.output,
    tf_fabric=tf_fabric,
    annotation_obj_specs=annotation_obj_specs,
    label_specs=label_specs,
    label_processors=processors,
)

print('RUNNING AUTOLABELER...')
labels = labeler.labelize()

# output labels to annotation sheet
sheet = BasicAnnotationSheet(
    annotations=labels,
    tf_fabric=tf_fabric,
)
print('BUILDING ANNOTATION SHEET...')
sheet.save_docx(str(snakemake.output))
