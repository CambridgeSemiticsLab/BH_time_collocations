"""This module contains a class for assigning labels to a set of Text-Fabric nodes."""

from abc import ABC, abstractmethod
from typing import Set, List, Dict
from tf.fabric import Fabric

from labeling.specifiers import LingLabel, ValueQuery


class BaseLabeler(ABC):
    """Base object for label processing."""

    @abstractmethod
    def label(
            self,
            annotation_objects: Dict[str, Set[int]],
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
            targets: Set[int],
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
        # run all value queries and collect their output
        print('\tRunning feature value queries...')
        labeled_targets: List[LingLabel] = []
        for value_query in self.value_queries:
            for target in value_query.value.label.targets:

                # get the nodes on which to run the value query
                target_set = annotation_objects.get(target.name)
                if not target_set:
                    raise Exception(f'Target {target.name} not identified by any query!')

                # execute the query and process results
                query_results = self._run_query(value_query.query, target_set)
                label_tuple = (value_query.value.label.name, value_query.value.name)
                print(f'\t\t{label_tuple}, {len(query_results)} results')
                for node in query_results:
                    labeled_targets.append(
                        LingLabel(
                            label=value_query.value.label.name,
                            value=value_query.value.name,
                            node=node,
                            target=target.name,
                        )
                    )
        return labeled_targets
