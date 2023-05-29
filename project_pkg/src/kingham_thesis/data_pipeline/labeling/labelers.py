"""This module contains a class for assigning labels to a set of Text-Fabric nodes."""

import pickle

from abc import ABC, abstractmethod
from typing import Set, List, Dict, Tuple
from tf.fabric import Fabric

from kingham_thesis.data_pipeline.labeling.specifiers import (
    ValueQuery, NodeIdentifier, LingLabel
)


class BaseLabeler(ABC):
    """Base object for label processing."""

    @abstractmethod
    def label(
            self,
            annotation_objects: Dict[str, Set[int]],
            custom_sets: Dict[str, Set[int]],
    ) -> List[LingLabel]:
        """
        Process targets for a given label name.

        :param annotation_objects: a dictionary with object names
            as keys, and a set of node integers as values
        :return: a list of LingLabel named tuples
        """


class QueryLabeler(BaseLabeler):
    """Processor for autolabeling with Text-Fabric queries."""

    def __init__(
            self,
            tf_fabric: Fabric,
            value_queries: List[ValueQuery],
    ) -> None:
        """
        Initialize the labeler.

        :param tf_fabric: Text-Fabric object to use for running the queries
        :param value_queries: a dictionary that maps labels to label values
            to queries
        :return: None
        """
        self.tf_fabric = tf_fabric
        self.api = tf_fabric.api
        self.value_queries = value_queries

    def _run_query(
            self,
            query: str,
            sets: Dict[str, Set[int]],
    ) -> Set[int]:
        """Execute query and return results."""
        result_set = self.api.S.search(
            query,
            shallow=True,
            sets=sets,
        )
        return result_set

    def label(
            self,
            annotation_objects: Dict[str, Set[int]],
            custom_sets: Dict[str, Set[int]],
    ) -> List[LingLabel]:
        """Assign labels to targets based on queries."""
        # run all value queries and collect their output
        print('\tRunning feature value queries...')
        # map to dict in order to de-dup multiple matches on label ids (i.e. sans value),
        # this will give priority to the last-defined matching pattern
        labeled_targets: Dict[Tuple[str, NodeIdentifier, str], LingLabel] = {}
        for value_query in self.value_queries:
            for target in value_query.value.label.targets:

                # get the nodes on which to run the value query
                if target.name not in annotation_objects:
                    raise Exception(f'Target {target.name} not identified by any query!')

                # execute the query and process results
                query_results = self._run_query(
                    query=value_query.query,
                    sets={**annotation_objects, **custom_sets},
                )
                label_tuple = (value_query.value.label.name, value_query.value.name)
                print(f'\t\t{label_tuple}, {len(query_results)} results')
                for node in query_results:
                    found_label = LingLabel(
                        label=value_query.value.label.name,
                        value=value_query.value.name,
                        nid=NodeIdentifier.from_node(node, self.api),
                        target=target.name,
                    )
                    labeled_targets[found_label.id] = found_label
        found_labels = list(labeled_targets.values())
        return found_labels
