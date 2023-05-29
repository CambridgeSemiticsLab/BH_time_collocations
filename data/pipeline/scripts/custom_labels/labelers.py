"""This module contains a class for assigning labels to a set of Text-Fabric nodes."""

import pickle

from typing import Set, List, Dict
from tf.fabric import Fabric

from kingham_thesis.data_pipeline.labeling.specifiers import NodeIdentifier, LingLabel
from kingham_thesis.data_pipeline.labeling.labelers import BaseLabeler


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
            custom_sets: Dict[str, Set[int]],
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
