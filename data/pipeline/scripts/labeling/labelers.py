"""This module contains a class for assigning labels to a set of Text-Fabric nodes."""

import pickle

from abc import ABC, abstractmethod
from typing import Set, List, Dict
from tf.fabric import Fabric

from labeling.specifiers import ValueQuery, NodeIdentifier, LingLabel


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
    ) -> List[LingLabel]:
        """Assign labels to targets based on queries."""
        # run all value queries and collect their output
        print('\tRunning feature value queries...')
        labeled_targets: List[LingLabel] = []
        for value_query in self.value_queries:
            for target in value_query.value.label.targets:

                # get the nodes on which to run the value query
                if target.name not in annotation_objects:
                    raise Exception(f'Target {target.name} not identified by any query!')

                # execute the query and process results
                query_results = self._run_query(value_query.query, annotation_objects)
                label_tuple = (value_query.value.label.name, value_query.value.name)
                print(f'\t\t{label_tuple}, {len(query_results)} results')
                for node in query_results:
                    labeled_targets.append(
                        LingLabel(
                            label=value_query.value.label.name,
                            value=value_query.value.name,
                            nid=NodeIdentifier.from_node(node, self.api),
                            target=target.name,
                        )
                    )
        return labeled_targets


class EnglishTenseLabeler(BaseLabeler):
    """Class to load and apply tense tags based on English translations."""

    TARGET_NAME = 'verb'
    LABEL_NAME = 'tense'
    TENSE_KEY = 'esv_TAMsimp'
    TENSE_MAP = {
        'PAST': 'past',
        '?PAST': 'past',
        "PAST PERF": 'past perf',
        'PRES PERF': 'pres perf',
        'PAST PROG': 'past prog',
        'PRES': 'pres',
        '?PRES': 'pres',
        'PRES PROG': 'pres prog',
        'FUT': 'fut',
        'FUT PROG': 'fut prog',
        'HAB used to': 'hab',
    }

    def __init__(
            self,
            tf_fabric: Fabric,
            tense_file: str,
    ) -> None:
        """Initialize the labeler."""
        self.tf_fabric = tf_fabric
        self.tf_api = self.tf_fabric.api
        self.tense_data = self._load_tense_file(tense_file)

    @staticmethod
    def _log(message):
        """Log a message."""
        print(f'\t{message}')

    @staticmethod
    def _load_tense_file(path: str):
        """Load the tense data from disk."""
        with open(path, 'rb') as infile:
            return pickle.load(infile)

    def _map_tense_value(self, tense: str):
        """Map an english tense value to a dataset value."""
        if "MOD" in tense:
            return "mod"
        else:
            return self.TENSE_MAP.get(tense)

    def _get_tense_tag(self, verb: int):
        """Attempt to assign a tense tag to a verbnode."""
        if self.tf_api.F.vt.v(verb) == 'impv':
            return 'impv'
        tense_data = self.tense_data.get(verb)
        if tense_data:
            tense_tag = tense_data[self.TENSE_KEY]
            mapped_tag = self._map_tense_value(tense_tag)
            return mapped_tag

    def label(
            self,
            annotation_objects: Dict[str, Set[int]],
    ) -> List[LingLabel]:
        """Label verb objects."""
        self._log('Running verb-tense labeler...')

        # catch when no target defined
        if self.TARGET_NAME not in annotation_objects:
            print('No "verb" targets defined! Returning empty list.')
            return []

        labels = []
        for verb in annotation_objects[self.TARGET_NAME]:
            tense_tag = self._get_tense_tag(verb)
            if tense_tag:
                labels.append(
                    LingLabel(
                        label=self.LABEL_NAME,
                        value=tense_tag,
                        nid=NodeIdentifier.from_node(verb, self.tf_api),
                        target=self.TARGET_NAME,
                    )
                )
        self._log(f'\t{len(labels)} tenses autolabeled')
        return labels
