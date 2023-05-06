"""This module contains a class for applying labels using Text-Fabric queries."""

from textwrap import dedent
from typing import Set, List, Dict
from tf.fabric import Fabric

from labeling.utils import BaseLabelProcessor, LingLabel


class LabelQuery:
    """Class for storing label-value queries."""

    def __init__(
            self,
            targets: Set[str],
            label: str,
            value: str,
            query: str,
    ) -> None:
        """
        Initialize the label query object.

        :param targets: a set of target object names
            to run the queries against
        :param label: the label for which the query
            is to assign a value
        :param value: the value assigned by executing
            the query
        :param query: the query itself, with a special
            `target` object in the query, which corresponds
            with a node stored in a target set found from
            targets, as defined in the targets param
        """
        self.label = label
        self.value = value
        self.targets = targets
        self._query = query

    @property
    def query(self) -> str:
        """Retrieve query value."""
        return dedent(self._query)


class QueryLabeler(BaseLabelProcessor):
    """Processor for autolabeling with Text-Fabric queries."""

    def __init__(
            self,
            tf_fabric: Fabric,
            label_queries: List[LabelQuery],
    ) -> None:
        """
        Initialize the labeler.

        :param tf_fabric: Text-Fabric object to use for running the queries
        :param label_queries: a dictionary that maps labels to label values
            to queries
        :return: None
        """
        self.tf_fabric = tf_fabric
        self.api = tf_fabric.api
        self.label_queries = label_queries

    def _run_query(
            self,
            query: str,
            targets: Set[str],
    ) -> Set[int]:
        """Execute query and return results."""
        result_set = self.api.S.search(
            query,
            shallow=True,
            sets={'target': targets},
        )
        return result_set

    def label(
            self,
            annotation_objects: Dict[str, Set[int]],
    ) -> List[LingLabel]:
        """Assign labels to targets based on queries."""
        labeled_targets: List[LingLabel] = []
        for label_query in self.label_queries:
            for target in label_query.targets:
                target_set = annotation_objects[target]
                print(f'Running query for: {label_query.label}={label_query.value}...')
                query_results = self._run_query(label_query.query, target_set)
                print(f'\tresults: {len(query_results)}')
                for node in query_results:
                    labeled_targets.append(
                        LingLabel(
                            label_query.label,
                            label_query.value,
                            node,
                            target,
                        )
                    )
        return labeled_targets
