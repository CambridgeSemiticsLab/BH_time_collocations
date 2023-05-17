"""Module for configuring the autolabelling process for this project."""

import collections
import docx
import json

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Dict, Optional, Type
from tf.fabric import Fabric

from labeling.labelers import BaseLabeler, QueryLabeler
from labeling.specifiers import (
    TargetSpec, TargetQuerySpecifier, LabelSpec, ValueSpec, ValueQuery, LingLabel
)
from labeling import annotation_sheets
from labeling.annotation_sheets import (
    BaseAnnotationSheet, BasicAnnotationSheet, AnnotationSheetSpecs
)


def get_available_sheets() -> Dict[str, Type[BaseAnnotationSheet]]:
    """Go through annotation sheet module and collect possible sheet classes."""
    sheet_map: Dict[str, Type[BaseAnnotationSheet]] = {}
    for variable in dir(annotation_sheets):
        type_instance = getattr(annotation_sheets, variable)
        if hasattr(type_instance, 'NAME'):
            sheet_map[type_instance.NAME] = type_instance
    return sheet_map


SHEET_MAP = get_available_sheets()


class BaseLabelingProject(ABC):
    """A base label project object."""

    # Fill in configs in child classes
    NAME = 'base'
    CONFIGS = {}

    def __init__(
            self,
            annotation_outdir: str,
            annotation_indir: str,
            archive_dir: str,
            tf_fabric: Fabric,
            extra_labelers: Optional[List[BaseLabeler]],
    ):
        """Setup Text-Fabric variables."""
        self.annotation_outdir = Path(annotation_outdir)
        self.annotation_indir = Path(annotation_indir)
        self.archive_dir = archive_dir
        self.tf_fabric = tf_fabric
        self.extra_labelers = extra_labelers or []
        self.target_specs: Dict[str, TargetSpec] = {}
        self.label_specs: Dict[str, LabelSpec] = {}
        self.value_specs: Dict[str, ValueSpec] = {}
        self._load_configs()

    @property
    def name(self) -> str:
        """Define a name for the project."""
        return self.NAME

    @property
    def random_seed(self) -> int:
        """Get a random seed for consistent sampling."""
        return 42

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

    def format_annotation_filepath(self, id_int: int):
        """Return docx filename."""
        return self.annotation_outdir / f'{self.name}_{id_int}.docx'

    def write_annotation_sheets(self, labels: List[LingLabel]) -> None:
        """Write annotation sheet to disk, to outdir."""
        # group all labels by sheet
        sheet_grouped_labels: Dict[str, List[LingLabel]] = collections.defaultdict(list)
        for label in labels:
            sheet = self.label_specs[label.label].sheet
            sheet_grouped_labels[sheet].append(label)

        # write sheets
        for i, (sheet_name, labels) in enumerate(sheet_grouped_labels.items(), 1):
            sheet_class: Type[BaseAnnotationSheet] = SHEET_MAP[sheet_name]
            sheet = sheet_class(
                annotations=labels,
                tf_fabric=self.tf_fabric,
                project_name=self.name,
            )
            filepath = self.format_annotation_filepath(i)
            sheet.to_docx(filepath)
            print(f'\tannotation sheet written to {filepath}')

    def read_annotation_sheets(self) -> Dict[Path, BaseAnnotationSheet]:
        """Read a completed annotation sheet from indir."""
        sheets_to_collect = set(
            spec.sheet for spec in self.label_specs.values()
        )
        sheets = {}
        for sheet_path in sorted(self.annotation_indir.glob('*.docx')):
            doc = docx.Document(sheet_path)
            doc_meta: AnnotationSheetSpecs = json.loads(doc.core_properties.comments)
            should_get_sheet = (
                doc_meta['project'] == self.name
                and doc_meta['sheet'] in sheets_to_collect
            )
            if should_get_sheet:
                sheet_class: Type[BaseAnnotationSheet] = SHEET_MAP[doc_meta["sheet"]]
                sheets[sheet_path] = sheet_class.from_doc(
                        document=doc,
                        tf_fabric=self.tf_fabric,
                        project_name=self.name,
                    )
        return sheets

    def _load_configs(self) -> None:
        """Read configs file from YAML and populate configs."""
        if not self.CONFIGS:
            raise NotImplementedError("Must implement CONFIGS!")
        for target in self.CONFIGS['targets']:
            self.target_specs[target] = TargetSpec(name=target)
        for label, label_data in self.CONFIGS['labels'].items():
            self.label_specs[label] = label_spec = LabelSpec(
                name=label,
                targets=set(
                    TargetSpec(target)
                    for target in label_data['targets']
                ),
                value_strings=set(label_data['values']),
                sheet=label_data['sheet'],
            )
            for value in label_data['values']:
                self.value_specs[value] = ValueSpec(
                    name=value,
                    label=label_spec,
                )


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


class TestLabelingProject(BaseLabelingProject):
    """Define a thesis-level labeling project."""

    NAME = "test"
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
                    "clause_x",
                ],
                "sheet": BasicAnnotationSheet.NAME,
            },
            "aspect": {
                "targets": [
                    "time_clause",
                ],
                "values": [
                    "ach_di",
                    "acc_in",
                    "sta_in",
                ],
                "sheet": BasicAnnotationSheet.NAME,
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
                    "2.1",
                ],
                "sheet": BasicAnnotationSheet.NAME,
            },
            "tense": {
                "targets": [
                    "verb",
                ],
                "values": [
                    "past",
                ],
                "sheet": BasicAnnotationSheet.NAME,
            },
        },
    }

    @property
    def target_queries(self) -> List[TargetQuerySpecifier]:
        """Define queries for identifying target nodes."""
        return [
            TargetQuerySpecifier(
                self.target_specs["time_clause"],
                TIMECLAUSE_QUERY,
                0.1,
            ),
            TargetQuerySpecifier(
                self.target_specs["time_phrase"],
                TIMEPHRASE_QUERY,
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
                    t:target
                    /with/
                        phrase function=Time
                        < phrase function=Pred
                    /-/
                """,
            ),
            #  --- aspect ---
            ValueQuery(
                self.value_specs["acc_in"],
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
