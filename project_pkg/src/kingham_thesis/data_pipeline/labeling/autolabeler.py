"""Module for autolabeling linguistic objects."""

import collections
import random

from typing import List, Dict, Set, Any, Iterable
from datetime import datetime
from tf.fabric import Fabric

from kingham_thesis.data_pipeline.labeling.specifiers import (
    TargetQuerySpecifier, LabelSpec, NodeIdentifier, LingLabel
)
from kingham_thesis.data_pipeline.labeling.projects import (
    BaseLabelingProject, SetFinder
)


class AutoLabeler:
    """Object for assigning labels to linguistic objects automatically."""

    def __init__(
            self,
            tf_fabric: Fabric,
            project: BaseLabelingProject,
    ) -> None:
        """Initialialize the autolabeler."""
        self.tf_fabric = tf_fabric
        self.tf_api = tf_fabric.api
        self.project = project
        self.labels_to_do = self._get_labels_todo_by_target(project.label_specs.values())

    @staticmethod
    def _get_labels_todo_by_target(label_specs: Iterable[LabelSpec]) -> Dict[str, Set[str]]:
        """Return a mapping between a target string and all labels to-do."""
        target_to_labels = {}
        for label_spec in label_specs:
            for target in label_spec.targets:
                target_to_labels.setdefault(target.name, set()).add(label_spec.name)
        return target_to_labels

    @staticmethod
    def _log(message: Any, ts=False, indent=1):
        """Print log messages."""
        indent_str = '\t' * indent
        now = f'{datetime.now()}  ' if ts else ''
        print(f'{indent_str}{now}{message}')

    def _collect_custom_sets(self) -> Dict[str, Set[int]]:
        """Collect custom sets using SetFinder functions."""
        custom_sets = {}
        for set_finder in self.project.set_finders:
            custom_sets.update(set_finder(self.tf_api))
        return custom_sets

    def _run_object_query(
            self,
            query: str,
            target_sets: Dict[str, Set[int]],
            custom_sets: Dict[str, Set[int]],
    ) -> Set[int]:
        """Run a Text-Fabric query for an annotation object."""
        result_set = self.tf_api.S.search(
            query,
            shallow=True,
            sets={**target_sets, **custom_sets},
        )
        return result_set

    def _sample_query_results(
            self,
            results: Set[int],
            sample_size: int,
    ) -> Set[int]:
        """Get a fractional sample of a query result."""
        random.seed(self.project.random_seed)
        sample = random.sample(sorted(results), sample_size)
        return set(sample)

    def _collect_annotation_objects(
            self,
            object_specs: List[TargetQuerySpecifier],
            custom_sets: Dict[str, Set[int]],
    ) -> Dict[str, Set[int]]:
        """Collect all annotation objects."""
        target_objects: Dict[str, Set[int]] = collections.defaultdict(set)
        for spec in object_specs:
            self._log(f'Executing query for TARGET={spec.target.name}')
            query_results = self._run_object_query(
                query=spec.query,
                target_sets=target_objects,
                custom_sets=custom_sets,
            )
            self._log(f'raw query results: {len(query_results)}', indent=2)
            if not spec.sample:
                target_objects[spec.target.name].update(query_results)
            else:
                sample = self._sample_query_results(query_results, spec.sample)
                self._log(f'subsampled to: {len(sample)}', indent=2)
                target_objects[spec.target.name].update(sample)
        return target_objects

    def _get_auto_labels(
            self,
            annotation_objects: Dict[str, Set[int]],
            custom_sets: Dict[str, Set[int]],
    ) -> List[LingLabel]:
        """Get autolabels for all targeted nodes."""
        # collect all labels
        auto_labels: List[LingLabel] = []
        covered_targets: Dict[str, Set[NodeIdentifier]] = collections.defaultdict(set)

        # collect all labels produced by processors
        for labeler in self.project.labelers:
            labels = labeler.label(annotation_objects, custom_sets)
            for label in labels:
                covered_targets[label.label].add(label.nid)
                auto_labels.append(label)

        # append empty labels for unlabeled targets
        n_unlabeled = collections.Counter()
        for target, nodes in annotation_objects.items():
            expected_labels = self.labels_to_do[target]
            for node in nodes:
                nid = NodeIdentifier.from_node(node, self.tf_api)
                for label_str in expected_labels:
                    if nid not in covered_targets[label_str]:
                        n_unlabeled[label_str] += 1
                        auto_labels.append(
                            LingLabel(
                                label=label_str,
                                value='',
                                nid=nid,
                                target=target,
                            )
                        )

        # give report on labeling outcome
        self._log('Empty Labels')
        self._log(n_unlabeled.most_common(), ts=False, indent=2)

        # done
        return auto_labels

    def labelize(self) -> List[LingLabel]:
        """Generate labels and output an annotation file."""
        custom_sets = self._collect_custom_sets()
        annotation_objects = self._collect_annotation_objects(self.project.target_queries, custom_sets)
        auto_labels = self._get_auto_labels(annotation_objects, custom_sets)
        return auto_labels
