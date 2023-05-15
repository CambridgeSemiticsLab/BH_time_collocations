"""Module to handle archiving and updating of accepted labels."""

import os
import json
import collections

from typing import List, Set
from pathlib import Path
from tf.fabric import Fabric
from labeling.specifiers import LingLabel, ArchivableLingLabel, NodeIdentifier
from labeling.projects import BaseLabelingProject


class LabelArchivist:
    """Object to archive labels and keep them in-sync with latest defined targets."""

    JSON_INDENT = 2

    def __init__(
            self,
            tf_fabric: Fabric,
            project: BaseLabelingProject,
    ):
        """Initialize the LabelArchivist."""
        self.tf_fabric = tf_fabric
        self.tf_api = self.tf_fabric.api
        self.project = project

    @staticmethod
    def _log(message):
        """Print a log message."""
        print(f'\t{message}')

    @property
    def archive_filepath(self) -> Path:
        """Get a formatted archive filepath."""
        filename = f'{self.project.name}.json'
        filepath = Path(self.project.archive_dir) / filename
        return filepath

    def _convert_labels_to_archivable(
            self,
            ling_labels: List[LingLabel]
    ) -> List[ArchivableLingLabel]:
        """Convert a list of LingLabel objects to ArchivableLingLabel objects."""
        return [
            ArchivableLingLabel.from_ling_label(label, self.tf_api)
            for label in ling_labels
        ]

    def _read_archive_file(self) -> Set[ArchivableLingLabel]:
        """Read in Archived nodes."""
        with open(self.archive_filepath, 'r') as infile:
            labels = set(
                ArchivableLingLabel.from_serialization(label)
                for label in json.load(infile)
            )
            self._log(f'{len(labels)} labels read from {self.archive_filepath}')
            return labels

    @staticmethod
    def _label_is_filled(label: LingLabel):
        """Check whether a supplied label is filled in."""
        return bool(label.value)

    def label_is_well_formed(self, label: LingLabel) -> bool:
        """Check whether a label conforms to the project definitions."""
        return (
            label.label in self.project.label_specs
            and label.value in self.project.label_specs[label.label].value_strings
            and label.target in self.project.target_specs
        )

    def _read_annotation_sheets(self) -> List[ArchivableLingLabel]:
        """Read fresh archive in from annotation sheet."""
        annotation_sheets = self.project.read_annotation_sheets()
        annotations = []
        ill_formed = []
        for path, sheet in annotation_sheets.items():
            sheet_annotations = []
            for annotation in sheet.annotations:
                if self._label_is_filled(annotation):
                    if self.label_is_well_formed(annotation):
                        sheet_annotations.append(
                            ArchivableLingLabel.from_ling_label(
                                annotation, self.tf_api
                            )
                        )
                    else:
                        ill_formed.append(annotation)
            self._log(f'{len(sheet_annotations)} labels read from {path}')
            if ill_formed:
                self._log(f'\t!! {len(ill_formed)} ill-formed annotations were ignored !!')
                print('\n'.join('\t\t'+str(label) for label in ill_formed))
            annotations.extend(sheet_annotations)
        return annotations

    def _get_archive(self) -> Set[ArchivableLingLabel]:
        """Get archive labels."""
        archive = set()
        if self.archive_filepath.exists():
            archive.update(self._read_archive_file())
        archive.update(self._read_annotation_sheets())
        self._log(f'{len(archive)} unique labels from all sources')
        self._log('')
        return archive

    def _nid_has_changed(self, node_id: NodeIdentifier) -> bool:
        """
        Check to see if a NodeIdentifier still reflects the underlying corpus.

        Since we expect node numbers to be mutable across corpus versions,
        we depend on the NodeIdentifier, which is a two-tuple of (otype, oslots)
        for an annotated node. If the `oslots` of the represented node in the NID
        change, we know that the label needs to be redone.
        """
        if node_id.otype == 'word':
            # NB: slots should never change between corpus updates,
            # which the logic here depends on
            return False
        current_node = self.tf_api.L.u(
            node_id.oslots[0],  # use first slot as ref point for node lookup
            node_id.otype,
        )[0]
        current_nid = NodeIdentifier(
            otype=node_id.otype,
            oslots=self.tf_api.L.d(current_node, "word")
        )
        return current_nid != node_id

    def _sync_archive_with_corpus(
            self,
            labels: Set[ArchivableLingLabel]
    ) -> Set[ArchivableLingLabel]:
        """Sync archived objects with the corpus."""
        return set(
            label for label in labels
            if not self._nid_has_changed(label.nid)
        )

    def _sync_archive_with_latest(
            self,
            archive: Set[ArchivableLingLabel],
            latest_labels_archivable: List[ArchivableLingLabel],
    ) -> Set[ArchivableLingLabel]:
        """Remove any archived labels no longer in the raw project results."""
        latest_label_map = {
            label.id: label
            for label in latest_labels_archivable
        }
        good = set()
        obsolete = collections.Counter()
        for archived_label in archive:

            # decide whether archived label has bee obsoleted
            latest_label = latest_label_map.get(archived_label.id)
            if not latest_label:
                is_good = False
            elif self._label_is_filled(latest_label):
                is_good = (archived_label == latest_label)
            else:
                is_good = (archived_label.id == latest_label.id)

            # add to set depending on status
            if is_good:
                good.add(archived_label)
            else:
                obsolete[(archived_label.label, archived_label.value)] += 1

        if obsolete:
            self._log(f'\t!! OBSOLETE LABELS PRUNED FROM ARCHIVE: !!')
            self._log(f'\t{obsolete.most_common()}')

        return good

    @staticmethod
    def _get_labels_to_do(
            archive: Set[ArchivableLingLabel],
            latest: List[ArchivableLingLabel],
    ) -> List[ArchivableLingLabel]:
        """Get diff of latest and archive to get to-do labels."""
        archived_ids = set(
            label.id for label in archive
        )
        return [
            label for label in latest
            if label.id not in archived_ids
        ]

    def _write_archive(self, curated_labels: Set[ArchivableLingLabel]) -> None:
        """Write curated ling labels to the archive."""
        sorted_archive = sorted(curated_labels)
        with open(self.archive_filepath, 'w') as outfile:
            json.dump(sorted_archive, outfile, indent=self.JSON_INDENT)
        self._log(f'Archived {len(curated_labels)} labels to {self.archive_filepath}')

    def _convert_archivable_to_label(
            self,
            labels: List[ArchivableLingLabel],
    ) -> List[LingLabel]:
        """Convert a list of ArchivableLingLabel objects to LingLabel objects."""
        # return as LingLabel objects
        return [
            LingLabel.from_archivable_label(label, self.tf_api)
            for label in labels
        ]

    def _preserve_annotation_sheet(self):
        """Change annotation sheet status to complete and set it to read-only."""


    def curate_collection(
            self,
            latest_labels: List[LingLabel]
    ) -> List[LingLabel]:
        """Curate the existing collection of linguistic labels and return a todo list."""
        latest_labels_archivable = self._convert_labels_to_archivable(latest_labels)
        archive = self._get_archive()
        self._log(f'archive length before syncing: {len(archive)}')
        archive = self._sync_archive_with_corpus(archive)
        self._log(f'archive length after corpus-sync: {len(archive)}')
        archive = self._sync_archive_with_latest(archive, latest_labels_archivable)
        self._log(f'archive length after project-sync: {len(archive)}')
        self._log('')
        self._log(f'labels-to-do before archive-sync: {len(latest_labels_archivable)}')
        labels_to_do = self._get_labels_to_do(archive, latest_labels_archivable)
        self._log(f'labels-to-do after archive-sync: {len(labels_to_do)}')
        self._log('')
        self._write_archive(archive)
        return self._convert_archivable_to_label(labels_to_do)
