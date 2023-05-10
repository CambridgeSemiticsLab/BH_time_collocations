"""Module for generating and ingesting annotation sheets."""

import collections
import docx

from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional
from docx.shared import Pt, RGBColor
from docx import Document
from docx.enum.text import WD_PARAGRAPH_ALIGNMENT
from docx.enum.style import WD_STYLE_TYPE
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from tf.fabric import Fabric
from labeling.utils import LingLabel


def add_hyperlink(paragraph, url, text, color=None, underline=True):
    """
    A function that places a hyperlink within a paragraph object.

    Source: https://github.com/python-openxml/python-docx/issues/74#issuecomment-261169410

    :param paragraph: The paragraph we are adding the hyperlink to.
    :param url: A string containing the required url
    :param text: The text displayed for the url
    :return: The hyperlink object
    """

    # This gets access to the document.xml.rels file and gets a new relation id value
    part = paragraph.part
    r_id = part.relate_to(url, docx.opc.constants.RELATIONSHIP_TYPE.HYPERLINK, is_external=True)

    # Create the w:hyperlink tag and add needed values
    hyperlink = docx.oxml.shared.OxmlElement('w:hyperlink')
    hyperlink.set(docx.oxml.shared.qn('r:id'), r_id, )

    # Create a w:r element
    new_run = docx.oxml.shared.OxmlElement('w:r')

    # Create a new w:rPr element
    rPr = docx.oxml.shared.OxmlElement('w:rPr')

    # Add color if it is given
    if not color is None:
        c = docx.oxml.shared.OxmlElement('w:color')
        c.set(docx.oxml.shared.qn('w:val'), color)
        rPr.append(c)

    # Remove underlining if it is requested
    if not underline:
        u = docx.oxml.shared.OxmlElement('w:u')
        u.set(docx.oxml.shared.qn('w:val'), 'none')
        rPr.append(u)

    # Join all the xml elements together add add the required text to the w:r element
    new_run.append(rPr)
    new_run.text = text
    hyperlink.append(new_run)

    paragraph._p.append(hyperlink)

    return hyperlink


class AnnotationSheet(ABC):
    """Object for generating an annotation sheet as a MS document."""

    def __init__(
            self,
            annotations: List[LingLabel],
            tf_fabric: Fabric,
    ):
        """Initialize an AnnotationSheet object."""
        self.annotations = annotations
        self.tf_fabric = tf_fabric
        self.tf_api = tf_fabric.api
        self.document = Document()
        self.styles = self._add_styles(self.document)
        self._build_document(self.document)

    @staticmethod
    @abstractmethod
    def _add_styles(self) -> Dict[str, Any]:
        """Set all styles for the document."""

    @abstractmethod
    def _build_document(self, document: Document()) -> None:
        """Build up document."""

    def save_docx(self, filepath: str) -> None:
        """Save a new docx with document."""
        self.document.save(filepath)

    @staticmethod
    @abstractmethod
    def _label_from_row(row) -> LingLabel:
        """Extract cell values from a table row into a LingLabel object."""

    @classmethod
    def from_docx(cls, filepath: str, tf_fabric: Fabric):
        """Read in docx annotation sheet."""
        document = docx.Document(filepath)
        annotations = []
        for table in document.tables:
            for row in table.rows:
                annotations.append(cls._label_from_row(row))
        return cls(annotations, tf_fabric)


class BasicAnnotationSheet(AnnotationSheet):
    """Prepare an annotation sheet."""

    GRAY = RGBColor(130, 130, 130)
    RED = RGBColor(255, 0, 0)
    SHEBANQ_LINK = (
        'https://shebanq.ancient-data.org/hebrew/text'
        '?book={book}&chapter={chapter}&verse={verse}&version=2021'
    )

    @staticmethod
    def _add_styles(document: Document()) -> Dict[str, Any]:
        """Set styles for the basic annotation sheet."""
        # dict to hold all new styles
        style_dict = {}
        styles = document.styles

        # set reference header style
        style_dict['ref'] = styles.add_style('Reference', WD_STYLE_TYPE.PARAGRAPH)
        style_dict['ref'].font.name = 'Times New Roman'
        style_dict['ref'].font.bold = True
        style_dict['ref'].font.size = Pt(12)
        style_dict['ref'].paragraph_format.alignment = WD_PARAGRAPH_ALIGNMENT.CENTER
        style_dict['ref'].paragraph_format.keep_with_next = True

        # set hebrew text style
        style_dict['hebrew'] = styles.add_style('Hebrew', WD_STYLE_TYPE.PARAGRAPH)
        style_dict['hebrew'].font.size = Pt(15)
        style_dict['hebrew'].font.name = 'SBL BibLit'
        style_dict['hebrew'].paragraph_format.alignment = WD_PARAGRAPH_ALIGNMENT.CENTER
        style_dict['hebrew'].paragraph_format.keep_with_next = True
        style_dict['hebrew'].paragraph_format.space_after = 0

        # done
        return style_dict

    def _get_clause_node(self, node: int) -> int:
        """Assign a clause node for a given node."""
        label_otype = self.tf_api.F.otype.v(node)
        clause_rank = self.tf_api.Nodes.otypeRank['clause']
        rank = self.tf_api.Nodes.otypeRank[label_otype]
        if rank > clause_rank:
            raise Exception(f'node {node} has a otype > clause!')
        elif label_otype == 'clause':
            return node
        else:
            return self.tf_api.L.u(node, 'clause')[0]

    def _cluster_labels_by_clause(self, labels: List[LingLabel]):
        """Cluster labels by clause."""
        cl_clustered_labels = collections.defaultdict(list)
        for label in labels:
            cl_node = self._get_clause_node(label.node)
            cl_clustered_labels[cl_node].append(label)
        return cl_clustered_labels

    def _add_reference_header(self, doc: Document, clause: int):
        """Add reference header to each entry."""
        book, ch, vs = self.tf_api.T.sectionFromNode(clause)
        ref = f'{clause}, {book} {ch}:{vs}'
        shebanq_link = self.SHEBANQ_LINK.format(
            book=book, chapter=str(ch), verse=str(vs)
        )
        heading = doc.add_paragraph(style=self.styles['ref'].name)
        add_hyperlink(heading, shebanq_link, ref)

    def _get_verse_nodes_from_clause(self, clause: int):
        """
        Get verse nodes from a clause node.

        Some clauses can span verses. So standard tf_api.L.u(clause, 'verse')
        will not work. We need to go from the slot level to get the one-to-many
        possible clause-verse relations.
        """
        verses = set()
        for word in self.tf_api.L.d(clause, 'word'):
            verses.add(self.tf_api.L.u(word, 'verse')[0])
        return sorted(verses)

    def _add_verse_text(self, doc: Document, clause: int):
        """Add verse hebrew text."""
        verse_nodes = self._get_verse_nodes_from_clause(clause)
        verse_text = doc.add_paragraph(
            self.tf_api.T.text(verse_nodes),
            style=self.styles['hebrew'].name,
        )
        verse_text.runs[0].font.color.rgb = self.GRAY

    def _add_clause_text(self, doc: Document, clause: int):
        """Add clause Hebrew text."""
        text = self.tf_api.T.text(clause)
        doc.add_paragraph(text, style=self.styles['hebrew'].name)

    def _add_slot_id_text(self, doc: Document, clause: int):
        """Add Hebrew text with slot identifiers indicated."""
        node_text = doc.add_paragraph(style=self.styles['hebrew'].name)
        for word in self.tf_api.L.d(clause, 'word'):
            node = node_text.add_run(str(word))
            node.font.superscript = True
            node.font.rtl = True
            node.font.color.rgb = RGBColor(255, 0, 0)
            heb_word = node_text.add_run(self.tf_api.T.text(word))
            heb_word.font.rtl = True
            heb_word.font.name = 'Times New Roman'
            heb_word.font.color.rgb = self.GRAY

    @staticmethod
    def _fix_autofit_bug(table):
        """
        Fix autofit bug for tables.

        Source:
        https://github.com/python-openxml/python-docx/issues/209#issuecomment-344417132
        """
        for col in table.columns:
            for cell in col.cells:
                cell._tc.tcPr.tcW.type = 'auto'

    def _get_label_node_text(self, label: LingLabel):
        """Get text for a label obj."""
        big_nodes = {
            'clause', 'sentence', 'verse', 'chapter', 'book'
        }
        otype = self.tf_api.F.otype.v(label.node)
        if otype in big_nodes:
            return f'[{otype}]'
        else:
            return self.tf_api.T.text(label.node)

    def _add_annotation_table(self, doc: Document, labels: List[LingLabel]):
        """Add annotation table to the document."""
        table = doc.add_table(rows=0, cols=5, style='Table Grid')
        table.alignment = WD_PARAGRAPH_ALIGNMENT.CENTER
        table.style.font.name = 'Helvetica Neue'
        table.style.paragraph_format.keep_with_next = True
        table.autofit = True
        for label in labels:
            row_cells = table.add_row().cells
            node_text = self._get_label_node_text(label)
            annotation_row = (
                label.label,
                str(label.node),
                label.target,
                node_text,
                label.value
            )
            for cell, text in zip(row_cells, annotation_row):
                cell.text = text
        self._fix_autofit_bug(table)

    def _build_document(self, document: Document()) -> None:
        """Build a document."""
        annotation_clusters = self._cluster_labels_by_clause(self.annotations)
        for clause in sorted(annotation_clusters):
            labels = annotation_clusters[clause]
            self._add_reference_header(document, clause)
            self._add_clause_text(document, clause)
            self._add_slot_id_text(document, clause)
            self._add_verse_text(document, clause)
            self._add_annotation_table(document, labels)
            document.add_paragraph('\n')

    @staticmethod
    def _label_from_row(row) -> LingLabel:
        """Extract a label from a row."""
        label_cell, node_cell, target_cell, text_cell, value_cell = row.cells
        ling_label = LingLabel(
            label=label_cell.text,
            value=value_cell.text,
            node=int(node_cell.text),
            target=target_cell.text
        )
        return ling_label
