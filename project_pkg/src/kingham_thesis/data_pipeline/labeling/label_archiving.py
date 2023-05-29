"""Module to handle archiving and updating of accepted labels."""

import json

from typing import List, Set, Dict
from pathlib import Path
from tf.fabric import Fabric
from kingham_thesis.data_pipeline.labeling.specifiers import LingLabel, NodeIdentifier
from kingham_thesis.data_pipeline.labeling.projects import BaseLabelingProject


class LabelArchivist:
    """Object to archive labels and keep them in-sync with latest defined targets."""

    JSON_INDENT = None

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

    @staticmethod
    def _label_is_filled(label: LingLabel):
        """Check whether a supplied label is filled in."""
        return bool(label.value)

    def _eval_label_health(self, label: LingLabel) -> Dict[str, bool]:
        """Check whether a label conforms to the project definitions."""
        value_check = (
            (label.value in self.project.label_specs[label.label].value_strings)
            if self.project.label_specs[label.label]
            else True
        )
        well_formed_eval = {
            'label in label specs': label.label in self.project.label_specs,
            'value in value specs': value_check,
            'target in target specs': label.target in self.project.target_specs
        }
        return well_formed_eval

    def _read_annotation_sheets(self) -> List[LingLabel]:
        """Read fresh archive in from annotation sheet."""
        annotation_sheets = self.project.read_annotation_sheets()
        annotations = []
        ill_formed = []
        for path, sheet in annotation_sheets.items():
            sheet_annotations = []
            for annotation in sheet.annotations:
                if self._label_is_filled(annotation):
                    label_health = self._eval_label_health(annotation)
                    if all(label_health.values()):
                        sheet_annotations.append(annotation)
                    else:
                        ill_formed.append((annotation, label_health))
            self._log(f'{len(sheet_annotations)} labels read from {path}')
            if ill_formed:
                self._log(f'\t!! {len(ill_formed)} ill-formed annotations were ignored !!')
                log_message = []
                for label, health_report in ill_formed:
                    log_message.append('\t\t' + str(label))
                    log_message.append('\t\t\t' + str(health_report))
                print('\n'.join(log_message))
            annotations.extend(sheet_annotations)
        return annotations

    def _get_archive(self) -> Set[LingLabel]:
        """Get archive labels."""
        archive = set()
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
            labels: Set[LingLabel]
    ) -> Set[LingLabel]:
        """Sync archived objects with the corpus."""
        return set(
            label for label in labels
            if not self._nid_has_changed(label.nid)
        )

    def _sync_archive_with_latest(
            self,
            archive: Set[LingLabel],
            latest_labels_archivable: List[LingLabel],
    ) -> Set[LingLabel]:
        """Remove any archived labels no longer in the raw project results."""
        latest_label_map = {
            label.id: label
            for label in latest_labels_archivable
        }
        good = set()
        obsolete = []

        for archived_label in archive:

            # decide whether archived label has been obsoleted
            latest_label = latest_label_map.get(archived_label.id)
            if not latest_label:
                is_good = False
                reason = 'not latest_label'
            else:
                is_good = (archived_label.id == latest_label.id)
                reason = 'archived_label.id != latest_label.id'

            # add to set depending on status
            if is_good:
                good.add(archived_label)
            else:
                obsolete.append((reason, archived_label, latest_label))

        if obsolete:
            log_message = []
            for reason, archived_label, latest_label in obsolete:
                log_message.append('\t\t' + str(archived_label))
                log_message.append(f'\t\t\tREASON: {reason}')
                log_message.append('\t\t\tLATEST: ' + str(latest_label))
            self._log(f'\t!! {len(obsolete)} OBSOLETE LABEL(S) PRUNED FROM ARCHIVE !!')
            print('\n'.join(log_message))

            self._log(f'')

        return good

    @staticmethod
    def _get_labels_to_do(
            archive: Set[LingLabel],
            latest: List[LingLabel],
    ) -> List[LingLabel]:
        """Get diff of latest and archive to get to-do labels."""
        archived_ids = set(
            label.id for label in archive
        )
        return [
            label for label in latest
            if label.id not in archived_ids
        ]

    def _write_archive(self, curated_labels: Set[LingLabel]) -> None:
        """Write curated ling labels to the archive."""
        sorted_archive = sorted(curated_labels)
        with open(self.archive_filepath, 'w') as outfile:
            json.dump(sorted_archive, outfile, indent=self.JSON_INDENT)
        self._log(f'Archived {len(curated_labels)} labels to {self.archive_filepath}')

    def _preserve_annotation_sheet(self):
        """Change annotation sheet status to complete and set it to read-only."""

    def curate_collection(
            self,
            latest_labels: List[LingLabel]
    ) -> List[LingLabel]:
        """Curate the existing collection of linguistic labels and return a todo list."""
        archive = self._get_archive()
        self._log(f'archive length before syncing: {len(archive)}')
        archive = self._sync_archive_with_corpus(archive)
        self._log(f'archive length after corpus-sync: {len(archive)}')
        archive = self._sync_archive_with_latest(archive, latest_labels)
        self._log(f'archive length after project-sync: {len(archive)}')
        self._log('')
        self._log(f'labels-to-do before archive-sync: {len(latest_labels)}')
        labels_to_do = self._get_labels_to_do(archive, latest_labels)
        self._log(f'labels-to-do after archive-sync: {len(labels_to_do)}')
        self._log('')
        self._write_archive(archive)
        return labels_to_do
