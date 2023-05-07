"""This module contains helper objects and definitions."""

from abc import ABC, abstractmethod
from typing import List, Dict, Set, NamedTuple, Tuple
from textwrap import dedent


class TargetObjectSpecifier(NamedTuple):
    """Object of interest to collect for annotation."""
    name: str
    raw_query: str  # string for Text-Fabric to query to identify the object

    @property
    def query(self) -> str:
        """Return dedented query string."""
        return dedent(self.raw_query)


class ObjectIdentifier(NamedTuple):
    """Slot-based identifier for a linguistic object."""
    slots: Tuple[int]
    otype: str

    def serialize(self) -> Tuple[Tuple[int], str]:
        """Serialize the identifier."""
        return self.slots, self.otype

    @classmethod
    def from_serialization(cls, serialization: Tuple[Tuple[int], str]) -> 'ObjectIdentifier':
        """Get an ObjectIdentifier from a serialized format."""
        return cls(*serialization)


class LingLabel(NamedTuple):
    """Object for storing linguistic labels."""
    label: str
    value: str
    node: int
    target: str


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
