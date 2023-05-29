"""This module contains helper objects and definitions."""

from typing import NamedTuple, Tuple, Set, Optional, Any

from textwrap import dedent
from tf.core.api import Api


def format_query(query: str) -> str:
    """Format a query string."""
    return dedent(query)


class TargetSpec(NamedTuple):
    """Specifier tuple for a target linguistic object."""

    name: str


class TargetQuerySpecifier(NamedTuple):
    """Object of interest to collect for annotation."""

    target: TargetSpec
    raw_query: str  # string for Text-Fabric to query to identify the object
    sample: Optional[int]

    @property
    def query(self) -> str:
        """Return query string."""
        return format_query(self.raw_query)


class LabelSpec(NamedTuple):
    """Object for storing label configurations."""

    name: str
    targets: Set[TargetSpec]
    value_strings: Optional[Set[str]]
    sheet: str


class ValueSpec(NamedTuple):
    """Specifier for values."""

    name: str
    label: LabelSpec


class ValueQuery(NamedTuple):
    """Class for automatically identifying a label value."""

    value: ValueSpec
    raw_query: str

    @property
    def query(self) -> str:
        """Retrieve query string."""
        return format_query(self.raw_query)


class NodeIdentifier(NamedTuple):
    """Class for representing a linguistic node without using the node number."""

    otype: str
    oslots: Tuple[int, ...]

    @classmethod
    def from_node(cls, node: int, tf_api: Api):
        """Get a node identifier tuple from a node."""
        otype = tf_api.F.otype.v(node)
        oslots = (
            tf_api.L.d(node, 'word')
            if otype != 'word'
            else (node,)
        )
        return cls(otype, oslots)


class LingLabel(NamedTuple):
    """
    Object for long-term storage of linguistic labels.

    Since node numbers might change after annotations are completed,
    we store the final annotation data under the slots and otype associated
    with the original node.
    """

    label: str
    value: str
    nid: NodeIdentifier
    target: str

    @property
    def id(self) -> Tuple[str, NodeIdentifier, str]:
        """Return a tuple to uniquely identify this label, without the filled value."""
        return self.label, self.nid, self.target


class SpecsDict(dict):
    """Class for holding specs for objects of interest."""

    def __getitem__(self, item: Any) -> Any:
        """Retrieve item from specs dict."""
        try:
            return super().__getitem__(item)
        except KeyError:
            raise Exception(f'Spec value "{item}" is undefined in configs!')
