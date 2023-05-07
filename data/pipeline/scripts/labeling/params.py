"""Module for configuring the autolabelling process for this project."""

from textwrap import dedent

from labeling.utils import TargetObjectSpecifier
from labeling.query_labeler import LabelQuery


annotation_obj_specs = [
    TargetObjectSpecifier(
        "time_clause",
        """
        clause
        /with/
            phrase function=Time
        /-/
        """,
    ),
    TargetObjectSpecifier(
        "time_phrase",
        """
        phrase function=Time
        """
    ),
    TargetObjectSpecifier(
        "verb",
        """
        w:word pdp=verb
        /with/
        clause
            phrase function=Time
            w
        /-/
        """
    ),
]


label_specs = {
    'time_clause': {'cl_type', 'aspect'},
    'time_phrase': {'tp_cluster'},
    'verb': {'tense'},
}


label_queries = [
    # clause types
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
    # aspect
    LabelQuery(
        targets={'time_clause'},
        label='aspect',
        value='ach_di',
        query="""
        t:target
        /with/
            word pdp=verb lex=BW>[|HLK[|CWB[|QRB[|>MR[
            phrase function=Cmpl
                word lex=>L pdp=prep
        /-/
        """
    )
]
