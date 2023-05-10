"""Module for autolabeling linguistic objects."""

import collections

from typing import List, Dict, Set, Any
from datetime import datetime
from tf.fabric import Fabric

from labeling.utils import TargetObjectSpecifier, LingLabel, BaseLabelProcessor


class AutoLabeler:
    """Object for assigning labels to linguistic objects automatically."""

    def __init__(
            self,
            outdir: str,
            tf_fabric: Fabric,
            annotation_obj_specs: List[TargetObjectSpecifier],
            label_specs: Dict[str, Set[str]],
            label_processors: List[BaseLabelProcessor],
    ) -> None:
        """Initialialize the autolabeler."""
        self.outdir = outdir
        self.annotation_obj_specs = annotation_obj_specs
        self.label_specs = label_specs
        self.label_processors = label_processors
        self.tf_fabric = tf_fabric
        self.tf_api = tf_fabric.api

    @staticmethod
    def _log(message: Any, ts=False, indent=0):
        """Print log messages."""
        indent_str = '\t' * indent
        now = f'{datetime.now()}  ' if ts else ''
        print(f'{indent_str}{now}{message}')

    def _run_object_query(self, query: str) -> Set[int]:
        """Run a Text-Fabric query for an annotation object."""
        result_set = self.tf_api.S.search(query, shallow=True)
        return result_set

    def _collect_annotation_objects(
            self,
            object_specs: List[TargetObjectSpecifier]
    ) -> Dict[str, Set[int]]:
        """Collect all annotation objects."""
        annotation_objects: Dict[str, Set[int]] = collections.defaultdict(set)
        for spec in object_specs:
            annotation_objects[spec.name] = self._run_object_query(spec.query)
        return annotation_objects

    def _filter_labeled_objects(
            self,
            annotation_objects: Dict[str, Set[int]],
    ) -> Dict[str, Set[int]]:
        """Filter out objects that have already been annotated."""
        # TODO
        return annotation_objects

    def _get_auto_labels(
            self,
            annotation_objects: Dict[str, Set[int]],
    ) -> List[LingLabel]:
        """Get autolabels for all targeted nodes."""
        # collect all labels
        auto_labels: List[LingLabel] = []
        covered_targets = collections.defaultdict(set)
        n_labeled = collections.Counter()

        # collect all labels produced by processors
        for processor in self.label_processors:
            for label in processor.label(annotation_objects):
                covered_targets[label.label].add(label.node)
                auto_labels.append(label)
                n_labeled[label.label] += 1

        # append empty labels for unlabeled targets
        n_unlabeled = collections.Counter()
        for name, nodes in annotation_objects.items():
            expected_labels = self.label_specs[name]
            for node in nodes:
                for label_str in expected_labels:
                    if node not in covered_targets[label_str]:
                        n_unlabeled[label_str] += 1
                        auto_labels.append(
                            LingLabel(label_str, '', node, name)
                        )

        # give report on labeling outcome
        self._log('**** Successfully Autolabeled ****')
        self._log(n_labeled.most_common(), ts=False, indent=1)
        self._log('********* Needs Labels ***********')
        self._log(n_unlabeled.most_common(), ts=False, indent=1)

        # done
        return auto_labels

    def labelize(self) -> List[LingLabel]:
        """Generate labels and output an annotation file."""
        annotation_objects = self._collect_annotation_objects(self.annotation_obj_specs)
        new_annotation_objs = self._filter_labeled_objects(annotation_objects)
        auto_labels = self._get_auto_labels(new_annotation_objs)
        return auto_labels
