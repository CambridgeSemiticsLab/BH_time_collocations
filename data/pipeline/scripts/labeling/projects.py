"""Module for configuring the autolabelling process for this project."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Dict, Optional, Type
from tf.fabric import Fabric
from labeling.labelers import BaseLabeler, QueryLabeler
from labeling.specifiers import (
    TargetSpec, TargetQuerySpecifier, LabelSpec, ValueSpec, ValueQuery, LingLabel
)
from labeling.annotation_sheets import BaseAnnotationSheet, BasicAnnotationSheet


# define some standard, inter-project templates
TIMEPHRASE_QUERY = "phrase function=Time"
TIMECLAUSE_QUERY = """
    clause
    /with/
        phrase function=Time
    /-/
"""
VERB_QUERY = """
    w:word pdp=verb
    /with/
    clause
        phrase function=Time
        w
    /-/
"""


class BaseLabelingProject(ABC):
    """A base label project object."""

    # Fill in configs in child classes
    CONFIGS = {}

    def __init__(
            self,
            annotation_outdir: str,
            annotation_indir: str,
            tf_fabric: Fabric,
            extra_labelers: Optional[List[BaseLabeler]],
    ):
        """Setup Text-Fabric variables."""
        self.annotation_outdir = annotation_outdir
        self.annotation_indir = annotation_indir
        self.tf_fabric = tf_fabric
        self.extra_labelers = extra_labelers or []
        self.targets: Dict[str, TargetSpec] = {}
        self.labels: Dict[str, LabelSpec] = {}
        self.values: Dict[str, ValueSpec] = {}
        self._load_configs()

    @property
    @abstractmethod
    def name(self) -> str:
        """Return a project name."""

    @property
    @abstractmethod
    def annotation_sheet(self) -> Type[BaseAnnotationSheet]:
        """Return AnnotationSheets."""

    @property
    @abstractmethod
    def target_queries(self) -> List[TargetQuerySpecifier]:
        """Return target queries."""

    @property
    def label_value_queries(self) -> List[ValueQuery]:
        """Return label value queries (optional if using querying for autolabeling)."""
        return []

    @property
    def labelers(self) -> List[BaseLabeler]:
        """Return default QueryLabeler + optional labelers."""
        return (
            [QueryLabeler(self.tf_fabric, self.label_value_queries)]
            + self.extra_labelers
        )

    @property
    def annotation_file(self):
        """Return docx filename."""
        return f'{self.name}.docx'

    def write_annotation_sheet(self, labels: List[LingLabel]) -> None:
        """Return instanced annotation sheets."""
        filepath = Path(self.annotation_outdir) / self.annotation_file
        sheet = self.annotation_sheet(
            annotations=labels,
            tf_fabric=self.tf_fabric,
            project_name=self.name,
        )
        sheet.to_docx(filepath)

    def read_annotation_sheet(self) -> BaseAnnotationSheet:
        """Read an annotation sheet in from disk."""
        filepath = Path(self.annotation_indir) / self.annotation_file
        if filepath.exists():
            # return populated sheet
            return self.annotation_sheet.from_docx(
                filepath=filepath,
                tf_fabric=self.tf_fabric,
                project_name=self.name,
            )
        else:
            # return blank sheet
            return self.annotation_sheet(
                annotations=[],
                tf_fabric=self.tf_fabric,
                project_name=self.name,
            )

    def _load_configs(self) -> None:
        """Read configs file from YAML and populate configs."""
        if not self.CONFIGS:
            raise NotImplementedError("Must implement CONFIGS!")
        for target in self.CONFIGS['targets']:
            self.targets[target] = TargetSpec(name=target)
        for label, label_data in self.CONFIGS['labels'].items():
            self.labels[label] = label_spec = LabelSpec(
                name=label,
                targets=set(
                    TargetSpec(target)
                    for target in label_data['targets']
                )
            )
            for value in label_data['values']:
                self.values[value] = ValueSpec(
                    name=value,
                    label=label_spec,
                )


class TestLabelingProject(BaseLabelingProject):
    """Define a thesis-level labeling project."""

    NAME = "test_1.0"
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
                "values": [
                    "x_clause",
                ],
            },
            "aspect": {
                "targets": [
                    "time_clause",
                ],
                "values": [
                    "ach_di",
                ],
            },
            "tp_cluster": {
                "targets": [
                    "time_phrase",
                ],
                "values": [
                    "1.1.1.1",
                    "1.1.1.2.1",
                    "1.1.1.2.2",
                    "1.1.1.3",
                    "1.1.2",
                ],
            },
            "tense": {
                "targets": [
                    "verb",
                ],
                "values": [
                    "past",
                ],
            },
        },
    }

    @property
    def name(self) -> str:
        """Define a name for the project."""
        return self.NAME

    @property
    def annotation_sheet(self) -> Type[BaseAnnotationSheet]:
        """Return annotation sheet class."""
        return BasicAnnotationSheet

    @property
    def target_queries(self) -> List[TargetQuerySpecifier]:
        """Define queries for identifying target nodes."""
        return [
            TargetQuerySpecifier(
                self.targets["time_clause"],
                TIMECLAUSE_QUERY,
            ),
            TargetQuerySpecifier(
                self.targets["time_phrase"],
                TIMEPHRASE_QUERY,

            ),
            TargetQuerySpecifier(
                self.targets["verb"],
                VERB_QUERY,
            ),
        ]

    @property
    def label_value_queries(self) -> List[ValueQuery]:
        """Define label value queries."""
        return [
            # --- clause types ---
            ValueQuery(
                self.values["x_clause"],
                """
                    t:target
                    /with/
                        phrase function=Time
                        < phrase function=Pred
                    /-/
                """,
            ),
            #  --- aspect ---
            ValueQuery(
                self.values["ach_di"],
                """
                t:target
                /with/
                    word pdp=verb lex=BW>[|HLK[|CWB[|QRB[|>MR[
                    phrase function=Cmpl
                        word lex=>L pdp=prep
                /-/
                """
            )
        ]
