import re
from pathlib import Path
import matplotlib.pyplot as plt

class Exporter:
    """Provide export capabilities for notebooks."""

    def __init__(self, outdir, analysis_name=''):
        self.outdir = Path(outdir).joinpath(analysis_name)
        if not self.outdir.exists():
            self.outdir.mkdir()
        
    def get_filename(self, name, suffix=''):
        """Provide a formatted filename."""
        filename = self.outdir.joinpath(name)
        filename = Path(f'{filename}.{suffix}')
        return filename

    def latex_input(self, name, graphname, **kwargs):
        """Format latex code for inputting a graph."""
        kgs = {
            'pos': 'htbp!',
            'align': '\centering',
            'label': name,
            'width': '0.6',
            'height': '',
        }
        kgs.update(kwargs)
        latex = (
            f'\\begin{{figure}}[{kgs["pos"]}]\n'
            f'{kgs["align"]}\n'
            '\caption{'+kgs.get('caption', '')+'}\n'
            '\label{fig:'+kgs['label']+'}\n'
            '\includegraphics[\n'
            '    width='+kgs['width']+'\\textwidth,\n'
            '    height='+kgs['height']+'\\textheight,\n'
            '    keepaspectratio,\n'
            ']{'+graphname+'}\n'
            '\end{figure}'
        )
        return latex

    def plot(self, name, save_latex=False, texkwargs={}, **savekwargs):
        """Exports a plot to the output directory.

        On best-practices for Latex figure formats:
        https://tex.stackexchange.com/questions/136087/selecting-best-file-extension-for-graphics-figures-pictures
        """
        savekwargs2 = dict(
            bbox_inches='tight',
            format='pdf',
        )
        savekwargs2.update(**savekwargs)
        filename = self.get_filename(name, savekwargs2['format'])
        plt.savefig(filename, **savekwargs2)

        if save_latex:
            latex_code = self.latex_input(
                name,
                filename.name,
                **texkwargs,
            ) 
            latexfile = Path(self.get_filename(name, 'tex'))
            latexfile.write_text(latex_code)
            

    def table(self, df, name, adjustbox=False, **kwargs):
        """Exports a table to Latex format."""
        tabkwargs = dict(
            escape=True,
            label=f'table:{name}',
            position='htbp!',
        )
        tabkwargs.update(kwargs)
        fname = self.get_filename(name, 'tex')
        table = df.to_latex(**tabkwargs)

        # un-escape macro brackets;
        # while it is possible to disable escaping,
        # this causes problems when unpredictable text enters in
        table = re.sub(
            r'\\textbackslash (.*?)\\{(.*?)\\}', 
            r'\\\g<1>{\g<2>}', 
            table,
        )

        # insert an adjustbox command to keep 
        # wider tables on the page
        if adjustbox:
            table = re.sub(
                r'(\n\\begin{tabular}.*)', 
                r'\n\\begin{adjustbox}{width=1\\textwidth}\g<1>',
                table,
            )
            table = re.sub(
                r'(\n\\end{tabular}.*)',
                r'\g<1>\n\\end{adjustbox}',
                table,
            )
        
        # done.
        # export the table
        Path(fname).write_text(table)

    def text(self, string, name):
        """Export a simple string (or int) for insertion in latex doc."""
        fname = self.get_filename(name, 'tex')
        string = str(string) + '%' # prevent extra space from being added
        fname.write_text(string)

    def number(self, number, name, commafy=True, roundto=None):
        """Export a number value."""
        fname = self.get_filename(name, 'tex')
        number = round(number, roundto)
        if commafy:
            number = "{:,}".format(number) 
        self.text(number, name)
