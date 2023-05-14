"""Module to handle archiving and updating of accepted labels."""

from typing import List
from tf.fabric import Fabric
from labeling.specifiers import LingLabel
from labeling.projects import BaseLabelingProject


class LabelArchivist:
    """Object to archive labels and keep them in-sync with latest defined targets."""

    def __init__(
            self,
            tf_fabric: Fabric,
            project: BaseLabelingProject
    ):
        """Initialize the LabelArchivist."""
        self.tf_fabric = tf_fabric
        self.tf_api = self.tf_fabric.api
        self.project = project

    def curate_collection(
            self,
            ling_labels: List[LingLabel]
    ) -> List[LingLabel]:
        """Curate the existing collection of linguistic labels and return a todo list."""
        # TODO
        sheet = self.project.read_annotation_sheet()
        return ling_labels
