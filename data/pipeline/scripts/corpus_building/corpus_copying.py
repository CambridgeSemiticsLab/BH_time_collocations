"""Module for copying a text-fabric corpus."""

from copy import deepcopy
from tf.fabric import Fabric
from typing import Dict, Any


def _copy_meta_dicts(feature_dict):
    """Extract metakwargs."""
    return {
        feat: deepcopy(feat_obj.meta)
        for feat, feat_obj in feature_dict.items()
    }


def _copy_feature_dicts(feature_dict):
    """Extract feature dicts."""
    return {
        feat: dict(feat_obj.items())
        for feat, feat_obj in feature_dict.items()
    }


def get_copy_of_corpus(tf_fabric: Fabric) -> Dict[str, Any]:
    """Get a copy of a corpus's resources."""
    tf_api = tf_fabric.api
    node_features = _copy_feature_dicts(tf_api.F.__dict__)
    edge_features = _copy_feature_dicts(tf_api.E.__dict__)
    metadata = dict(
        **_copy_meta_dicts(tf_api.F.__dict__),
        **_copy_meta_dicts(tf_api.E.__dict__),
    )
    metadata['otext'] = tf_fabric.features['otext'].metaData
    return {
        'nodeFeatures': node_features,
        'edgeFeatures': edge_features,
        'metaData': metadata,
    }
