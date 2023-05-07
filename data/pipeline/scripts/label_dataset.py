"""This module labels the analysis dataset."""

from tf.fabric import Fabric

from labeling.autolabeler import AutoLabeler
from labeling.query_labeler import QueryLabeler
from labeling.params import (
    annotation_obj_specs, label_specs, label_queries,
)

tf_fabric = Fabric(
    locations=snakemake.input.corpus,
    silent="deep",
)
tf_api = tf_fabric.loadAll()


processors = [
    QueryLabeler(tf_fabric, label_queries),
]


labeler = AutoLabeler(
    outdir=snakemake.output,
    tf_fabric=tf_fabric,
    annotation_obj_specs=annotation_obj_specs,
    label_specs=label_specs,
    label_processors=processors,
)

labeler.labelize()
