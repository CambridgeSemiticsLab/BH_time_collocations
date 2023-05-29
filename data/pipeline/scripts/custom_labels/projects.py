"""Module for configuring the autolabelling process for this project."""


from typing import List, Dict, Type


from kingham_thesis.data_pipeline.labeling.specifiers import TargetQuerySpecifier, ValueQuery
from kingham_thesis.data_pipeline.labeling.projects import BaseLabelingProject
from kingham_thesis.data_pipeline.labeling.annotation_sheets import BaseAnnotationSheet

from custom_labels.annotation_sheets import BasicAnnotationSheet


# define some standard, inter-project templates
TIMECLAUSE_QUERY = """
    time_clause:clause
    /with/
        phrase function=Time
    /-/
"""
TIMEPHRASE_QUERY = """
    time_phrase:phrase function=Time
    /with/
    time_clause
        time_phrase
    /-/
"""
VERB_QUERY = """
    verb:word pdp=verb
    /with/
    time_clause
        phrase function=Time
        verb
    /-/
"""


CL_TYPE_VALUES = [
    "x_clause",
    "clause_x",
    "x_clause_x",
    "wayehi_x",
    "wehaya_x",
    "medial",
    "nmcl",
    "nmcl_disloc",  # nominal clause dislocation, e.g. Num 10:10
    "ellp",         # ellipsis
]

ASPECT_VALUES = [
    "sta_tr",
    "sta_ac",
    "sta_in",
    "sta_po",
    "ach_rd",
    "ach_id",
    "ach_cy",
    "act_di",
    "act_un",
    "acc_in",
    "acc_ru",
    "acc_di_iter",
    "none",  # none, not applicable
]

TP_CLUSTER_VALUES = [
    "1.1.1.2.1.1",
    "1.1.1.2.1.1.1",
    "1.1.1.2.1.2",
    "1.1.1.2.1.2.1",
    "1.1.1.2.2",
    "1.1.1.2.5",
    "1.1.1.2.2.1",
    "1.1.1.1",
    "1.1.2.3",
    "1.1.1.3",
    "1.1.2.4",
    "1.1.2.1.3",
    "1.1.2.2",
    "1.1.2.1.1",
    "1.1.1.2.1.1.2",
    "1.1.2.5.1",
    "1.1.2.5.2",
    "1.1.2.5.3.1",
    "1.1.2.5.3.1.1",
    "1.1.2.5.3.2",
    "1.1.2.6",
    "1.1.2.7",
    "1.1.1.2.1.2.2",
    "1.1.1.2.1.1.3",
    "1.1.1.2.2.2",
    "1.1.1.2.4.1",
    "1.1.1.2.4.1.1",
    "1.1.2.4.1",
    "1.1.2.5.3.2.1",
    "1.1.2.5.3.2.2",
    "1.1.2.5.3.2.3",
    "1.1.2.5.3.2.1.1",
    "1.1.2.5.3.2.1.2",
    "1.1.2.6.1",
    "1.1.1.2.3",
    "1.1.1.2.2.3",
    "1.1.2.1.1.1",
    "1.1.2.1.3.1",
    "1.1.2.1.3.2",
    "1.1.1.2.2.4",
    "1.1.2.5.3.2.1.3",
    "1.1.2.1.2",
]

TENSE_VALUES = [
    "past",  # preterite
    "pres perf",
    "past perf",
    "past prog",
    "pres",
    "pres prog",
    "fut",
    "fut prog",
    "mod",
    "epis mod",
    'impv',
    "gnom",  # gnomic
    "hab",  # iterative / habitual
    "inf",  # infinitival
    "ptcp",  # participial
]


class BTimeLabelingProject(BaseLabelingProject):
    """Define a thesis-level labeling project."""

    NAME = "b_time"
    CONFIGS = {
        "targets": [
            "time_clause",
            "time_phrase",
            "verb",
        ],
        "labels": {
            "cl_type": {
                "targets": [
                    "time_clause",
                ],
                "values": CL_TYPE_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
            "aspect": {
                "targets": [
                    "time_clause",
                ],
                "values": ASPECT_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
            "tp_cluster": {
                "targets": [
                    "time_phrase",
                ],
                "values": TP_CLUSTER_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
            "tense": {
                "targets": [
                    "verb",
                ],
                "values": TENSE_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
        },
    }
    LABEL_ORDER = {
        "cl_type": 0,
        "aspect": 1,
        "tense": 2,
        "tp_cluster": 3,
        "tp_head": 4,
    }

    @property
    def sheet_map(self) -> Dict[str, Type[BaseAnnotationSheet]]:
        """Return mapping from sheet name to its constructor class."""
        return {
            BasicAnnotationSheet.NAME: BasicAnnotationSheet
        }

    @property
    def target_queries(self) -> List[TargetQuerySpecifier]:
        """Define queries for identifying target nodes."""
        return [
            TargetQuerySpecifier(
                self.target_specs["time_phrase"],
                """
                phrase function=Time
                /with/
                    =: word lex=B
                /or/
                    =: word pdp=advb
                    <: word lex=B
                /-/
                """,
                None,
            ),
            TargetQuerySpecifier(
                self.target_specs["time_clause"],
                """
                clause
                /with/
                    time_phrase
                /-/
                """,
                None,
            ),
            TargetQuerySpecifier(
                self.target_specs["verb"],
                VERB_QUERY,
                None,
            ),
        ]

    @property
    def label_value_queries(self) -> List[ValueQuery]:
        """Define label value queries."""
        return [
            # --- clause types ---
            ValueQuery(
                self.value_specs["x_clause"],
                """
                    time_clause
                    /with/
                        phrase function=Time
                        < phrase function=Pred|PreO|PreS
                    /-/
                """,
            ),
            ValueQuery(
                self.value_specs["clause_x"],
                """
                    time_clause
                    /without/
                        time_phrase
                        < phrase function=Objc|Cmpl|Subj|Pred|PreO|PreS
                    /-/
                    /with/
                        phrase function=Pred|PreO|PreS
                            w:word pdp=verb vt=wayq|perf|impf
                            /without/
                            w lex=HJH[ vt=wayq nu=sg ps=p3
                            /-/
                            /without/
                            word lex=W
                            <: w lex=HJH[ vt=perf nu=sg ps=p3
                            /-/
                        < time_phrase
                    /-/
                """
            ),
            ValueQuery(
                self.value_specs["medial"],
                """
                time_clause
                /with/
                    phrase function=Pred|PreO|PreS
                    < time_phrase
                    < phrase function=Subj|Objc|Cmpl
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["wayehi_x"],
                """
                time_clause
                /with/
                    word lex=HJH[ vt=wayq nu=sg ps=p3
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["wehaya_x"],
                """
                time_clause
                /with/
                    word lex=W
                    <: word lex=HJH[ vt=perf nu=sg ps=p3
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["nmcl"],
                """
                time_clause typ=NmCl|AjCl
                """
            ),
            #  --- aspect ---
            ValueQuery(
                self.value_specs["acc_in"],
                """
                time_clause
                /with/
                    word pdp=verb lex=BW>[|HLK[|CWB[|QRB[|>MR[
                    phrase function=Cmpl
                        word lex=>L pdp=prep
                /-/
                """
            ),

            # --- TP Cluster --
            ValueQuery(
                self.value_specs["1.1.1.2.1.1"],  # distal demonstrative
                """
                time_phrase
                /with/
                    =: word lex=B
                    <: word lex=H
                /-/
                /with/
                    word lex=HJ>|HMH|HM|HW>
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.1.2.1.2"],  # proximal demonstrative
                """
                time_phrase
                /with/
                    =: word lex=B
                    <: word lex=H
                /-/
                /with/
                    word lex=Z>T|>LH|ZH
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.1.2.2"],  # ordinal
                """
                time_phrase
                /with/
                    =: word lex=B
                    <: word lex=H
                /-/
                /with/
                    word lex=H 
                    word ls=ordn
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.1.1"],  # definite, standalone
                """
                phrase function=Time
                /without/
                    word lex=HJ>|HMH|HM|HW>|Z>T|>LH|ZH
                /-/
                /without/
                    word ls=ordn|card
                /-/
                    =: word lex=B
                    <: w:word lex=H
                    /without/
                    phrase
                        w
                        < word lex=H
                    /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.2.1.1"],  # calendrical
                """
                time_phrase
                /with/
                    =: word lex=B
                    <: word ls=card
                    < word lex=L
                    <: word lex=H
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.2.5.1"],  # clause-anchored
                """
                tp:time_phrase
                /with/
                clause
                    tp
                        =: word lex=B
                        <: lastword:word st=c
                <: clause
                    
                tp := lastword
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.2.5.2"],  # construct, cardinal anchored
                """
                tp:time_phrase
                /with/
                    =: word lex=B
                    <: cons_word:word st=c ls#card lex#<WD/
                    <: word ls=card
                    last_word:word

                tp := last_word
                last_word # cons_word
                /-/
                """
            ),
            ValueQuery(
                self.value_specs["1.1.2.6"],  # suffix-anchored
                """
                tp:time_phrase
                /with/
                    =: word lex=B
                    <: last_word:word st=a prs#absent

                tp := last_word
                /-/
                """
            )
        ]


class NonTemporalClauseProject(BaseLabelingProject):
    """Project to label non-temporal clauses."""

    NAME = "nt_clauses"
    CONFIGS = {
        "targets": [
           "nt_clause",
           "verb",
        ],
        "labels": {
            "cl_type": {
                "targets": ["nt_clause"],
                "values": CL_TYPE_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
            "aspect": {
                "targets": ["nt_clause"],
                "values": ASPECT_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
            "tense": {
                "targets": ["verb"],
                "values": TENSE_VALUES,
                "sheet": BasicAnnotationSheet.NAME,
            },
        },
    }
    LABEL_ORDER = {
        "cl_type": 0,
        "aspect": 1,
        "tense": 2,
    }

    @property
    def target_queries(self) -> List[TargetQuerySpecifier]:
        """Define queries for identifying target nodes."""
        return [
            TargetQuerySpecifier(
                self.target_specs["nt_clause"],
                """
                """,
                1100,
            )
        ]
