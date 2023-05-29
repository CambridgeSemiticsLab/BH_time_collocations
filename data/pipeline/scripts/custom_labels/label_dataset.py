"""
A module for applying hand-checked labels to the dataset.

NB: This scripting module is to be executed by Snakemake, which
injects a variable, `snakemake`, for accessing program parameters.
"""

from tf.fabric import Fabric

from kingham_thesis.data_pipeline.labeling.project_runner import ProjectRunner

from projects import BTimeLabelingProject, NonTemporalClauseProject
from labelers import EnglishTenseLabeler


# configure resources
tf_fabric = Fabric(
    locations=snakemake.input.corpus,
    silent="deep",
)
tf_api = tf_fabric.loadAll()

# set up projects
base_project_params = dict(
    annotation_dir=snakemake.params.annotation_dir,
    tf_fabric=tf_fabric,
)

english_tense_labeler = EnglishTenseLabeler(
    tf_fabric, str(snakemake.input.tense_data)
)

projects = [
    BTimeLabelingProject(
        **base_project_params,
        extra_labelers=[english_tense_labeler],
    ),
    NonTemporalClauseProject(
        **base_project_params,
        extra_labelers=[english_tense_labeler],
    ),
]

projects_todo = [
    project for project in projects
    if project.name in snakemake.params.projects
]

# run projects
runner = ProjectRunner(projects_todo, tf_fabric)
runner.build_labels()
