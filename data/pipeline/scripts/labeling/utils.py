"""This module contains helper objects and definitions."""

from abc import ABC, abstractmethod
from typing import List, Dict, Set, NamedTuple, Tuple
from textwrap import dedent
from tf.core.api import Api


class TargetObjectSpecifier(NamedTuple):
    """Object of interest to collect for annotation."""

    name: str
    raw_query: str  # string for Text-Fabric to query to identify the object

    @property
    def query(self) -> str:
        """Return dedented query string."""
        return dedent(self.raw_query)


class LingLabel(NamedTuple):
    """Object for storing linguistic labels."""

    label: str
    value: str
    node: int
    target: str


class NodeIdentifier(NamedTuple):
    """Class for representing a linguistic node without using the node number."""

    otype: str
    oslots: Tuple[int]


class FrozenLingLabel(NamedTuple):
    """
    Object for long-term storage of linguistic labels.

    Since node numbers might change after annotations are completed,
    we store the final annotation data under the slots and otype associated
    with the original node. Obsoleted FrozenLingLabels can easily be
    identified by searching for nodes with the same otype and oslots
    within a newer version of the corpus. Failure to find a match indicates
    the label should be redone.
    """

    label: str
    value: str
    nid: NodeIdentifier
    target: str

    @classmethod
    def from_ling_label(cls, ling_label: LingLabel, tf_api: Api) -> 'FrozenLingLabel':
        """Get FrozenLingLabel from LingLabel object."""
        otype = tf_api.F.otype.v(ling_label.node)
        oslots = (
            tf_api.L.d(ling_label.node, 'word')
            if otype != 'word'
            else ling_label.node
        )
        return FrozenLingLabel(
            label=ling_label.label,
            value=ling_label.value,
            nid=NodeIdentifier(otype, oslots),
            target=ling_label.target,
        )


class BaseLabelProcessor(ABC):
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
