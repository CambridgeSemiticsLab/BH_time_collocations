"""Module for configuring the autolabelling process for this project."""

from textwrap import dedent

from labeling.utils import TargetObjectSpecifier
from labeling.query_labeler import LabelQuery


annotation_obj_specs = [
    TargetObjectSpecifier(
        "time_clause",
        dedent("""
            clause
            /with/
                phrase function=Time
            /-/
        """),
    ),
    TargetObjectSpecifier(
        "time_phrase",
        dedent("""
            phrase function=Time
        """)
    ),
    TargetObjectSpecifier(
        "verb",
        dedent("""
            w:word pdp=verb
            /with/
            clause
                phrase function=Time
                w
            /-/
        """)
    ),
]


label_specs = {
    'time_clause': {'cl_type'},
    'time_phrase': {'tp_cluster'},
    'verb': {'tense'},
}


label_queries = [
    LabelQuery(
        targets={'time_clause'},
        label='cl_type',
        value='x_clause',
        query="""
            t:target
            /with/
                phrase function=Time
                < phrase function=Pred
            /-/
        """,
    ),
]
