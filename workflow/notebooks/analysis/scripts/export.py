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
        filename = f'{filename}.{suffix}'
        return filename

    def plot(self, name, **kwargs):
        """Exports a plot to the output directory.

        On best-practices for Latex figure formats:
        https://tex.stackexchange.com/questions/136087/selecting-best-file-extension-for-graphics-figures-pictures
        """
        save_kwargs = dict(
            bbox_inches='tight',
            format='pdf',
        )
        save_kwargs.update(**kwargs)
        filename = self.get_filename(name, save_kwargs['format'])
        plt.savefig(filename, **save_kwargs)

    def table(self, df, name, adjustbox=False, **kwargs):
        """Exports a table to Latex format."""
        tabkwargs = dict(
            escape=False,
            position='htbp!',
        )
        tabkwargs.update(kwargs)
        fname = self.get_filename(name, 'tex')
        table = df.to_latex(**tabkwargs)

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
