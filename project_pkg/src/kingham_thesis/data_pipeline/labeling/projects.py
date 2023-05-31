"""Module for configuring the autolabelling process for this project."""

import collections
import docx
import json

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Dict, Optional, Type, Any, Callable, Set
from tf.fabric import Fabric
from tf.core.api import Api

from kingham_thesis.data_pipeline.labeling.labelers import BaseLabeler, QueryLabeler
from kingham_thesis.data_pipeline.labeling.specifiers import (
    TargetSpec, TargetQuerySpecifier, LabelSpec,
    ValueSpec, ValueQuery, LingLabel, SpecsDict
)
from kingham_thesis.data_pipeline.labeling.annotation_sheets import (
    BaseAnnotationSheet, AnnotationSheetSpecs
)


SetFinder = Callable[[Api], Dict[str, Set[int]]]


class BaseLabelingProject(ABC):
    """A base label project object."""

    # Fill in configs in child classes
    NAME = 'base'
    CONFIGS = {}

    def __init__(
            self,
            annotation_dir: str,
            tf_fabric: Fabric,
            extra_labelers: Optional[List[BaseLabeler]] = None,
    ) -> None:
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
    @abstractmethod
    def sheet_map(self) -> Dict[str, Type[BaseAnnotationSheet]]:
        """Return a mapping from sheet names to their class types."""

    @property
    def random_seed(self) -> int:
        """Get a random seed for consistent sampling."""
        return 42

    @property
    def set_finders(self) -> List[SetFinder]:
        """Get a list of methods for finding set definitions."""
        return []

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

    def _format_annotation_filepath(self, id_int: int) -> Path:
        """Return docx filename."""
        return self.annotation_outdir / f'{self.name}_{id_int}.docx'

    def _write_annotation_metadata(self, metadata, id_int: int) -> None:
        """Write annotation metadata to json."""
        filepath = self.annotation_outdir / f'{self.name}_{id_int}_metadata.json'
        with open(filepath, 'w') as outfile:
            json.dump(metadata, outfile)

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
        # initialize the output directory
        self._initialize_outdir()

        # if no remaining labels, do nothing
        if not labels:
            print('\tNo more labels to review! No sheet written.')
            return

        # group all labels by sheet
        sheet_grouped_labels: Dict[str, List[LingLabel]] = collections.defaultdict(list)
        for label in labels:
            sheet = self.label_specs[label.label].sheet
            sheet_grouped_labels[sheet].append(label)

        # write sheets
        annotation_id_start = self._get_annotation_sheet_id_start()
        for i, (sheet_name, labels) in enumerate(sheet_grouped_labels.items(), annotation_id_start):
            sheet_class: Type[BaseAnnotationSheet] = self.sheet_map[sheet_name]
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
                sheet_class: Type[BaseAnnotationSheet] = self.sheet_map[doc_meta["sheet"]]
                sheets[sheet_path] = sheet_class.from_doc(
                    document=doc,
                    tf_fabric=self.tf_fabric,
                    project=self,
                    metadata=metadata,
                )
        return sheets

    def _load_configs(self) -> None:
        """Read configs from the class attributes and populate configs dicts."""
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
