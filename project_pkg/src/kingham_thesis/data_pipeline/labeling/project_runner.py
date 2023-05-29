"""Module contains classes for running projects."""

from typing import List
from tf.fabric import Fabric
from kingham_thesis.data_pipeline.labeling.projects import BaseLabelingProject
from kingham_thesis.data_pipeline.labeling.autolabeler import AutoLabeler
from kingham_thesis.data_pipeline.labeling.label_archiving import LabelArchivist


class ProjectRunner:
    """Object for running labeling projects."""

    def __init__(
            self,
            projects: List[BaseLabelingProject],
            tf_fabric: Fabric,
    ):
        """Initialize the ProjectRunner."""
        self.projects = projects
        self.tf_fabric = tf_fabric

    @staticmethod
    def _log(*message):
        """Basic logger with newline."""
        print()
        print(*message)

    def build_labels(self):
        """Run the project labels."""
        for project in self.projects:
            self._log(f"--- BEGIN {project.name} ---")

            self._log('RUNNING AUTOLABELER...')
            labeler = AutoLabeler(
                tf_fabric=self.tf_fabric,
                project=project,
            )
            labels = labeler.labelize()

            self._log('ARCHIVING COMPLETED LABELS & GETTING DIFFS...')
            archivist = LabelArchivist(
                tf_fabric=self.tf_fabric,
                project=project,
            )
            labels_to_do = archivist.curate_collection(labels)

            self._log('BUILDING ANNOTATION SHEET...')
            project.write_annotation_sheets(labels_to_do)

            self._log('--- DONE ---')
            self._log()
