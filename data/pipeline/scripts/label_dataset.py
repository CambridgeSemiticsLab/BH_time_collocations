"""
A module for applying hand-checked labels to the dataset.

NB: This scripting module is to be executed by Snakemake, which
injects a variable, `snakemake`, for accessing program parameters.
"""

from tf.fabric import Fabric
from labeling.project_runner import ProjectRunner
from labeling.projects import TestLabelingProject


# configure resources
tf_fabric = Fabric(
    locations=snakemake.input.corpus,
    silent="deep",
)
tf_api = tf_fabric.loadAll()

# set up projects
projects = [
    TestLabelingProject(
        annotation_outdir=str(snakemake.params.annotation_outdir),
        annotation_indir=str(snakemake.params.annotation_indir),
        archive_dir=str(snakemake.params.archive_dir),
        tf_fabric=tf_fabric,
        extra_labelers=[],
    ),
]
projects_todo = [
    project for project in projects
    if project.name in snakemake.params.projects
]


# run projects
runner = ProjectRunner(projects_todo, tf_fabric)
runner.build_labels()
