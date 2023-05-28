"""Module for configuring the autolabelling process for this project."""

import collections
import docx
import json

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Dict, Optional, Type, Any
from tf.fabric import Fabric

from labeling.labelers import BaseLabeler, QueryLabeler
from labeling.specifiers import (
    TargetSpec, TargetQuerySpecifier, LabelSpec,
    ValueSpec, ValueQuery, LingLabel, SpecsDict
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
            annotation_dir: str,
            tf_fabric: Fabric,
            extra_labelers: Optional[List[BaseLabeler]],
    ):
        """Setup Text-Fabric variables."""
        self.annotation_dir = Path(annotation_dir)
        self.sheets_dir = self.annotation_dir / "sheets"
        self.annotation_outdir = (
                self.sheets_dir / "blank" / self.name
        )
        self.annotation_indir = (
                self.sheets_dir / "complete" / self.name
        )
        self.archive_dir = self.annotation_dir / "json"
        self.tf_fabric = tf_fabric
        self.extra_labelers = extra_labelers or []
        self.target_specs: Dict[str, TargetSpec] = SpecsDict()
        self.label_specs: Dict[str, LabelSpec] = SpecsDict()
        self.value_specs: Dict[str, ValueSpec] = SpecsDict()
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

    def _format_annotation_filepath(self, id_int: int):
        """Return docx filename."""
        return self.annotation_outdir / f'{self.name}_{id_int}.docx'

    def _write_annotation_metadata(self, metadata, id_int: int) -> None:
        """Write annotation metadata to json."""
        filepath = self.annotation_outdir / f'{self.name}_{id_int}_metadata.json'
        with open(filepath, 'w') as outfile:
            json.dump(metadata, outfile, indent=2)

    def _read_annotation_metadata(self, filestem: str) -> Dict[str, Any]:
        """Read annotation metadata."""
        filepath = self.annotation_indir / f'{filestem}_metadata.json'
        with open(filepath, 'r') as infile:
            return json.load(infile)

    def _initialize_outdir(self) -> None:
        """Make the outdir if it doesn't exist."""
        if not self.annotation_outdir.exists():
            self.annotation_outdir.mkdir()

    def _get_annotation_sheet_id_start(self) -> int:
        """Retrieve the annotation sheet id based on completed sheets."""
        n_completed_sheets = len(list(
            self.annotation_indir.glob('*.docx')
        ))
        return n_completed_sheets + 1

    def write_annotation_sheets(self, labels: List[LingLabel]) -> None:
        """Write annotation sheet to disk, to outdir."""
        # group all labels by sheet
        sheet_grouped_labels: Dict[str, List[LingLabel]] = collections.defaultdict(list)
        for label in labels:
            sheet = self.label_specs[label.label].sheet
            sheet_grouped_labels[sheet].append(label)

        # write sheets
        self._initialize_outdir()
        annotation_id_start = self._get_annotation_sheet_id_start()
        for i, (sheet_name, labels) in enumerate(sheet_grouped_labels.items(), annotation_id_start):
            sheet_class: Type[BaseAnnotationSheet] = SHEET_MAP[sheet_name]
            sheet = sheet_class(
                annotations=labels,
                tf_fabric=self.tf_fabric,
                project=self,
            )
            filepath = self._format_annotation_filepath(i)
            sheet.to_docx(filepath)
            self._write_annotation_metadata(sheet.annotation_metadata, i)
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
                metadata = self._read_annotation_metadata(sheet_path.stem)
                sheet_class: Type[BaseAnnotationSheet] = SHEET_MAP[doc_meta["sheet"]]
                sheets[sheet_path] = sheet_class.from_doc(
                    document=doc,
                    tf_fabric=self.tf_fabric,
                    project=self,
                    metadata=metadata,
                )
        return sheets

    def _load_configs(self) -> None:
        """Read configs file from YAML and populate configs."""
        if not self.CONFIGS:
            raise NotImplementedError("Must implement CONFIGS!")
        for target in self.CONFIGS['targets']:
            self.target_specs[target] = TargetSpec(name=target)
        for label, label_data in self.CONFIGS['labels'].items():
            values = label_data.get('values')
            value_strings = set(values) if values else None
            self.label_specs[label] = label_spec = LabelSpec(
                name=label,
                targets=set(
                    TargetSpec(target)
                    for target in label_data['targets']
                ),
                value_strings=value_strings,
                sheet=label_data['sheet'],
            )
            for value in label_data.get('values', []):
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
                "values": [
                    "x_clause",
                    "clause_x",
                    "x_clause_x",
                    "wayehi_x",
                    "wehaya_x",
                    "medial",
                    "nmcl",
                    "nmcl_disloc",  # nominal clause dislocation, e.g. Num 10:10
                    "ellp",  # ellipsis
                ],
                "sheet": BasicAnnotationSheet.NAME,
            },
            "aspect": {
                "targets": [
                    "time_clause",
                ],
                "values": [
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
                ],
                "sheet": BasicAnnotationSheet.NAME,
            },
            "tp_cluster": {
                "targets": [
                    "time_phrase",
                ],
                "values": [
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
                ],
                "sheet": BasicAnnotationSheet.NAME,
            },
            "tense": {
                "targets": [
                    "verb",
                ],
                "values": [
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
                ],
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
