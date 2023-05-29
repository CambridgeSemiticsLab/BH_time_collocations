"""Provide classes and methods for constructing a corpus."""

import json
import shutil
import collections

from pprint import pformat
from textwrap import indent
from pathlib import Path
from typing import Set, Dict, List, Optional, Any, Tuple, Counter

from tf.fabric import Fabric

from kingham_thesis.data_pipeline.corpus_building.canonical_sorting import canonical_order
from kingham_thesis.data_pipeline.corpus_building.corpus_copying import get_copy_of_corpus
from kingham_thesis.data_pipeline.corpus_building.utils import (
    FeatureGenerator, MetaDataDict, NodeFeatureDict, EdgeFeatureDict,
    EditAction, CorpusData,
)


class ThesisCorpusBuilder:
    """Class for building thesis corpus."""

    def __init__(
            self,
            locations: Optional[List[str]] = None,
            filter_metadata: Optional[Dict[str, Any]] = None,
            annotation_files: Optional[List[Path]] = None,
            book_limit: Optional[int] = None,
            delete_features: Optional[Set[str]] = None,
            rename_features: Optional[Dict[str, str]] = None,
            add_features: Optional[Dict[str, FeatureGenerator]] = None,
            update_metadata: Optional[MetaDataDict] = None,
            update_features: Optional[NodeFeatureDict] = None,
            update_edges: Optional[EdgeFeatureDict] = None,
            delete_nodes: Optional[Set[int]] = None,
            edit_actions: Optional[List[EditAction]] = None,
            tf_fabric: Optional[Fabric] = None,
    ):
        """Initialize the thesis corpus builder."""
        self.locations = locations or ''
        self.filter_metadata = filter_metadata or {}
        self.annotation_files = annotation_files
        self.book_limit = book_limit
        self.add_features = add_features or {}
        self.update_metadata = update_metadata or {}
        self.delete_features = delete_features or set()
        self.rename_features = rename_features or {}
        self.delete_nodes = delete_nodes or set()
        self.update_features = update_features or {}
        self.update_edges = update_edges or {}
        self._add_edit_actions(
            edit_actions or [],
            self.delete_nodes,
            self.update_features,
            self.update_edges,
        )
        self.tf_fabric = tf_fabric
        self.tf_api = tf_fabric.api if tf_fabric else None

    @staticmethod
    def _add_edit_actions(
            edit_actions: List[EditAction],
            delete_nodes: Set[int],
            update_features: NodeFeatureDict,
            update_edges: EdgeFeatureDict,
    ) -> None:
        """Add all edit actions to the correct dicts / sets."""
        for action in edit_actions:
            delete_nodes.update(action.deletions)
            for feature_name, feature_dict in action.feature_updates.items():
                update_features.setdefault(feature_name, {}).update(feature_dict)
            for edge_name, edge_dict in action.edge_updates.items():
                update_edges.setdefault(edge_name, {}).update(edge_dict)

    def _get_keep_node_set(self):
        """Get set of nodes to keep."""
        keep_nodes = set()
        max_slot = 0
        book_limit = (
            self.tf_api.T.nodeFromSection((self.book_limit,))
            if self.book_limit else None
        )
        for book_node in self.tf_api.F.otype.s('book'):
            if book_limit and book_node > book_limit:
                break
            keep_nodes.add(book_node)
            for node in self.tf_api.L.d(book_node):
                keep_nodes.add(node)
                if self.tf_api.F.otype.v(node) == 'word':
                    max_slot = node
        return keep_nodes, max_slot

    @staticmethod
    def _filter_nodes_from_feature_dict(feature_dict, keep_nodes):
        """Filter keep nodes."""
        filtered_feature_dict = {}
        for feature, node_dict in feature_dict.items():
            filtered_feature_dict[feature] = {
                node: feature
                for node, feature in node_dict.items()
                if node in keep_nodes
            }
        return filtered_feature_dict

    def _filter_feature_dict_nodes(self, corpus_data, keep_node_set):
        """Filter feature dict nodes."""
        corpus_data['nodeFeatures'] = self._filter_nodes_from_feature_dict(
            corpus_data['nodeFeatures'],
            keep_node_set,
        )
        corpus_data['edgeFeatures'] = self._filter_nodes_from_feature_dict(
            corpus_data['edgeFeatures'],
            keep_node_set
        )

    def _rebuild_nodes_from_oslots(
            self,
            oslot_map,
            otype_map,
            max_slot,
    ) -> Dict[int, int]:
        """Rebuild node numbering scheme from oslots."""
        # get sorted list of oslot source_data
        oslots = []
        for node, oslot_set in oslot_map.items():
            otype = otype_map[node]
            otype_rank = self.tf_api.Nodes.otypeRank[otype]
            oslots.append((node, otype_rank, set(oslot_set)))
        oslots.sort(key=canonical_order)

        # create mapping to new node numbers
        new_node_map = {
            old_node: (i + max_slot)
            for i, (old_node, _, _) in enumerate(oslots, 1)
        }
        return new_node_map

    @staticmethod
    def _reindex_node_features(old_node_features, remapper):
        """Reindex node features."""
        node_features = collections.defaultdict(dict)
        for feature, node_dict in old_node_features.items():
            for node, fvalue in node_dict.items():
                node_features[feature][remapper(node)] = fvalue
        return node_features

    @staticmethod
    def _reindex_edge_features(old_edge_features, remapper):
        """Reindex edge features."""
        edge_features = collections.defaultdict(dict)
        for feature, edge_dict in old_edge_features.items():
            for node, edges in edge_dict.items():
                if isinstance(edges, dict):
                    edge_features[feature][remapper(node)] = {
                        remapper(n) for n, v
                        in edges.items()
                    }
                else:
                    edge_features[feature][remapper(node)] = set(
                        remapper(n) for n in edges
                    )
        return edge_features

    def _reindex_nodes(self, corpus_data: CorpusData, max_slot):
        """Reindex nodes."""
        # rebuild node numbering from oslot source_data
        new_node_map = self._rebuild_nodes_from_oslots(
            corpus_data['edgeFeatures']['oslots'],
            corpus_data['nodeFeatures']['otype'],
            max_slot,
        )

        # remap all node features using the new numbering scheme
        remapper = lambda node: new_node_map.get(node, node)
        node_features = self._reindex_node_features(corpus_data['nodeFeatures'], remapper)
        edge_features = self._reindex_edge_features(corpus_data['edgeFeatures'], remapper)

        # apply the changes to the dict
        corpus_data['nodeFeatures'] = node_features
        corpus_data['edgeFeatures'] = edge_features

        # add a helper map back to old nodes for referencing
        corpus_data['edgeFeatures']['omap@2021-KT'] = {
            new_node: {old_node: None}
            for old_node, new_node in new_node_map.items()
        }

    def _update_metadata(self, corpus_data: CorpusData):
        """Update metadata fields."""
        for field, data in self.update_metadata.items():
            corpus_data['metaData'].setdefault(field, {}).update(data)

    def _filter_metadata(self, corpus_data: CorpusData):
        """Remap metadata for this project."""
        # delete old generic metadata from all metadata dicts
        for feature, meta in corpus_data['metaData'].items():
            corpus_data['metaData'][feature] = {
                k: v for k, v in meta.items()
                if k not in self.filter_metadata
            }

    def _delete_features(self, corpus_data: CorpusData):
        """Remove features from the dataset."""
        for feature in self.delete_features:
            for data_type, data_dict in corpus_data.items():
                if feature in data_dict:
                    del data_dict[feature]

    def _rename_features(self, corpus_data: CorpusData):
        """Rename features in the dataset."""
        for old_name, new_name in self.rename_features.items():
            for data_type, data_dict in corpus_data.items():
                if old_name in data_dict:
                    data_dict[new_name] = data_dict[old_name]
                    del data_dict[old_name]

    def _update_features(self, corpus_data: CorpusData):
        """Update feature node mappings."""
        for feature, update_dict in self.update_features.items():
            corpus_data['nodeFeatures'].setdefault(feature, {}).update(update_dict)

    def _update_edges(self, corpus_data: CorpusData):
        """Update edge feature node mappings."""
        for feature, update_dict in self.update_edges.items():
            corpus_data['edgeFeatures'].setdefault(feature, {}).update(update_dict)

    def _delete_nodes_from_features(self, corpus_data: CorpusData):
        """Delete nodes from feature dicts."""
        for feature, node_data in corpus_data['nodeFeatures'].items():
            corpus_data['nodeFeatures'][feature] = {
                node: value
                for node, value in node_data.items()
                if node not in self.delete_nodes
            }

    def _delete_nodes_from_edges(self, corpus_data: CorpusData) -> None:
        """Delete nodes from edge relations."""
        # delete from edge relations
        new_edges = collections.defaultdict(dict)
        for feature, edge_data in corpus_data['edgeFeatures'].items():
            for node, edges in edge_data.items():
                if node in self.delete_nodes:
                    continue
                elif isinstance(edges, dict):
                    remaining_edges = {
                        n: value
                        for n, value in edges.items()
                        if n not in self.delete_nodes
                    }
                else:
                    remaining_edges = set(
                        n for n in edges
                        if n not in self.delete_nodes
                    )
                if remaining_edges:
                    # only add non-empty edge relations
                    new_edges[feature][node] = remaining_edges
        corpus_data['edgeFeatures'] = new_edges

    def _delete_nodes(self, corpus_data: CorpusData) -> None:
        """Delete nodes from the corpus."""
        self._delete_nodes_from_features(corpus_data)
        self._delete_nodes_from_edges(corpus_data)

    @staticmethod
    def _clear_directory(dest_dir: str):
        """Empty a destination directory of old source_data."""
        dest_dir = Path(dest_dir)
        if dest_dir.exists():
            shutil.rmtree(dest_dir)
        dest_dir.mkdir(parents=True)

    def _load_tf_corpus(self):
        """Load Text Fabric corpus."""
        if not self.tf_fabric:
            self.tf_fabric = Fabric(self.locations)
            self.tf_api = self.tf_fabric.loadAll()

    @staticmethod
    def _load_json(filepath: Path):
        """Load json from a filepath."""
        return json.loads(filepath.read_text())

    def _read_annotation_files(self) -> Dict[str, List[Any]]:
        """Read annotation data from disk."""
        file2labels = {}
        for filepath in self.annotation_files:
            file2labels[filepath] = self._load_json(filepath)
        return file2labels

    def _get_nid_to_node(
            self,
            corpus_data: CorpusData
    ) -> Dict[Tuple[str, Tuple[int, ...]], int]:
        """Build reverse index of node id (otype, oslots) to new node number."""
        nid_to_node = {}
        for node, oslots in corpus_data['edgeFeatures']['oslots'].items():
            otype = corpus_data['nodeFeatures']['otype'][node]
            sorted_oslots = tuple(sorted(oslots))  # ensure oslots are sorted
            nid = (otype, sorted_oslots)
            nid_to_node[nid] = node
            # add slots too
            for slot in sorted_oslots:
                nid_to_node[('word', (slot,))] = slot
        return nid_to_node

    @staticmethod
    def _report_annotated_features(
            added_labels: Counter[str],
            unable_to_assign: List[Any],
    ) -> None:
        """Report on annotated features."""
        print(f'\t{sum(added_labels.values())} new node features assigned...')
        print(indent(pformat(added_labels.most_common()), '\t\t'))
        if unable_to_assign:
            print(f'\t{len(unable_to_assign)} labels could not be found in corpus!')
            for label in unable_to_assign:
                print(f'\t\t{label}')

    def _add_annotation_features(self, corpus_data: CorpusData):
        """Add features from custom annotations to new nodes."""
        file2labels = self._read_annotation_files()
        if file2labels:
            print('\tAdding new features from the following sources:')
        else:
            print('\tFound no pending annotation files. Skipping.')
        added_labels = collections.Counter()
        unable_to_assign = []
        nid_to_node = self._get_nid_to_node(corpus_data)
        for file, labels in file2labels.items():
            print('\t\t', file)
            for label_data in labels:
                label, value, raw_nid, target = label_data
                if raw_nid[0] == 'word':
                    # hand slots differently since they have no oslots entry
                    node = raw_nid[1][0]
                    corpus_data['nodeFeatures'][label][node] = value
                    corpus_data['nodeFeatures']['target'][node] = target
                    added_labels[label] += 1
                    added_labels['target'] += 1
                else:
                    nid = (raw_nid[0], tuple(raw_nid[1]))
                    node = nid_to_node.get(nid)
                    if node:
                        # add the label to the corpus
                        corpus_data['nodeFeatures'][label][node] = value
                        corpus_data['nodeFeatures']['target'][node] = target
                        added_labels[label] += 1
                        added_labels['target'] += 1
                    else:
                        unable_to_assign.append(label_data)
        self._report_annotated_features(added_labels, unable_to_assign)

    def build(self, dest_dir: str):
        """Build the corpus."""
        print('Loading BHSA corpus...')
        self._load_tf_corpus()

        print('Getting a copy of the corpus...')
        corpus_data = get_copy_of_corpus(self.tf_fabric)

        print('Filtering the nodes...')
        keep_node_set, max_slot = self._get_keep_node_set()
        self._filter_feature_dict_nodes(corpus_data, keep_node_set)

        print('Applying graph edits...')
        self._delete_nodes(corpus_data)
        self._update_features(corpus_data)
        self._update_edges(corpus_data)

        print('Reindexing the nodes...')
        self._reindex_nodes(corpus_data, max_slot)

        print('Rebuilding metadata...')
        self._filter_metadata(corpus_data)
        self._update_metadata(corpus_data)

        print('Refactoring features...')
        self._delete_features(corpus_data)
        self._rename_features(corpus_data)

        print('Adding new annotation features...')
        self._add_annotation_features(corpus_data)

        print('Saving BHSA-KT...')
        self._clear_directory(dest_dir)
        saver = Fabric(dest_dir)
        saver.save(**corpus_data)
