"""This module contains helper objects and definitions."""

from typing import NamedTuple, Tuple, Set, List, Any, Optional
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
    sample: Optional[float]

    @property
    def query(self) -> str:
        """Return query string."""
        return format_query(self.raw_query)


class LabelSpec(NamedTuple):
    """Object for storing label configurations."""

    name: str
    targets: Set[TargetSpec]
    value_strings: Set[str]
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


class LingLabel(NamedTuple):
    """Object for storing linguistic labels."""

    label: str
    value: str
    node: int
    target: str

    @classmethod
    def from_archivable_label(
            cls,
            archivable_label: 'ArchivableLingLabel',
            tf_api: Api,
    ) -> 'LingLabel':
        """Get LingLabel from an ArchivableLingLabel (see below)"""
        if archivable_label.nid.otype != 'word':
            node = tf_api.L.u(
                archivable_label.nid.oslots[0],
                archivable_label.nid.otype,
            )[0]
        else:
            node = archivable_label.nid.oslots[0]
        return cls(
            label=archivable_label.label,
            value=archivable_label.value,
            node=node,
            target=archivable_label.target,
        )


class NodeIdentifier(NamedTuple):
    """Class for representing a linguistic node without using the node number."""

    otype: str
    oslots: Tuple[int]


class ArchivableLingLabel(NamedTuple):
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

    @property
    def id(self) -> Tuple[str, NodeIdentifier, str]:
        """Return a tuple to unique identify this label, without the filled value."""
        return self.label, self.nid, self.target

    @classmethod
    def from_ling_label(cls, ling_label: LingLabel, tf_api: Api) -> 'ArchivableLingLabel':
        """Get FrozenLingLabel from LingLabel object."""
        otype = tf_api.F.otype.v(ling_label.node)
        oslots = (
            tf_api.L.d(ling_label.node, 'word')
            if otype != 'word'
            else (ling_label.node,)
        )
        return ArchivableLingLabel(
            label=ling_label.label,
            value=ling_label.value,
            nid=NodeIdentifier(otype, oslots),
            target=ling_label.target,
        )

    @classmethod
    def from_serialization(cls, serialization: List[Any]) -> 'ArchivableLingLabel':
        """Read in as serialization."""
        label, value, raw_nid, target = serialization
        return cls(
            label,
            value,
            NodeIdentifier(
                raw_nid[0],
                tuple(raw_nid[1]),
            ),
            target
        )
