"""Helper objects and definitions."""

from typing import (
    Dict, Optional, Union, Callable,
    TypedDict, Set
)


# Text-Fabric expected types
FeatureType = Union[str, int]
EdgeType = Union[Set[int], Dict[int, FeatureType]]
NodeFeatureDict = Dict[str, Dict[int, FeatureType]]
EdgeFeatureDict = Dict[str, Dict[int, EdgeType]]
MetaDataDict = Dict[str, Dict[str, str]]


class CorpusData(TypedDict):
    """Object to hold corpus source_data."""
    nodeFeatures: NodeFeatureDict
    edgeFeatures: EdgeFeatureDict
    metaData: MetaDataDict


FeatureGenerator = Callable[[CorpusData], NodeFeatureDict]


class EditAction:
    """Class for grouping related corpus edits."""

    def __init__(
            self,
            edge_updates: Optional[EdgeFeatureDict] = None,
            feature_updates: Optional[NodeFeatureDict] = None,
            deletions: Optional[Set[int]] = None,
            description: str = '',
    ):
        """Initialize edit action object."""
        self.deletions = deletions or set()
        self.feature_updates = feature_updates or {}
        self.edge_updates = edge_updates or {}
        self.description = description
