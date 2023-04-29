import re
from pathlib import Path
import matplotlib.pyplot as plt
from textwrap import dedent
from analysis_tools.src.analysis_tools.text.show import get_spread

def wrap_hebtext(string):
    """Wrap hebrew text with Latex polyglossia tag."""
    return '\\texthebrew{%s}'%string

class Exporter:
    """Provide export capabilities for notebooks."""

    def __init__(self, outdir, analysis_name=''):
        self.outdir = Path(outdir).joinpath(analysis_name)
        if not self.outdir.exists():
            self.outdir.mkdir()
        
    def get_subdir(self, name):
        """Get a subdirectory; make one if needed."""
        dir = self.outdir.joinpath(name)
        if not dir.exists():
            dir.mkdir()
        return dir
    
    def get_filename(self, name, subdir, suffix=''):
        """Provide a formatted filename."""
        filename = subdir.joinpath(name)
        filename = Path(f'{filename}.{suffix}')
        return filename
    
    def latex_input(self, name, graphname, **kwargs):
        """Format thesis code for inputting a graph."""
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
        subdir = self.get_subdir('plots')
        filename = self.get_filename(
            name, 
            subdir,
            suffix=savekwargs2['format']
        )
        plt.savefig(filename, **savekwargs2)

        if save_latex == True:
            latex_code = self.latex_input(
                name,
                filename.name,
                **texkwargs,
            ) 
            latexfile = self.get_filename(
                    name,
                    subdir,
                    suffix='tex',
            )
            latexfile.write_text(latex_code)

    def table(self, df, name, adjustbox=False, hebaxis=None, hebcols=[], **kwargs):
        """Exports a table to Latex format."""
        tabkwargs = dict(
            escape=True,
            label=f'table:{name}',
            position='htbp!',
        )
        tabkwargs.update(kwargs)
        fname = self.get_filename(
            name, 
            self.get_subdir('tables'),
            suffix='tex'
        )
        
        # update table with thesis formatting around columns containing hebrew
        if hebaxis == 0:
            df.index = df.index.to_series().apply(wrap_hebtext)
        elif hebaxis == 1:
            df.columns = df.columns.to_series().apply(wrap_hebtext)
        if hebcols:
            for col in hebcols:
                df[col] = df[col].apply(wrap_hebtext)
        
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

        return df

    def text(self, string, name):
        """Export a simple string (or int) for insertion in thesis doc."""
        fname = self.get_filename(
            name, 
            self.get_subdir('text'),
            suffix='tex'
        )
        string = str(string) + '%' # prevent extra space from being added
        fname.write_text(string)

        return string

    def number(self, number, name, commafy=True, roundto=None):
        """Export a number value."""
        number = round(number, roundto)
        if commafy:
            number = "{:,}".format(number) 
        self.text(number, name)
        return number

    
    def examples(self, df, name, textcols=['clause'], 
                 refcol='verse', joiner=' ', spread=0):
        
        """Copy Latex-formatted text examples to clipboard."""

        orig_shape = df.shape
        if spread > 0:
            spread_i = get_spread(df.index, spread)
            df = df.loc[spread_i]
            print(f'exporting {df.shape[0]} of {orig_shape[0]}...')
        
        # merge text (if more than 1 col included)
        texts = df[textcols].agg(joiner.join, 1)
        exs = []
        for i,node in enumerate(texts.index):      
            ref = df[refcol][node]
            text = dedent('''
            \\texthebrew{
            %s
            } (%s)
            ''') % (texts[node], ref)
            label = f'{name}{i+1}'
            exs.append('\\ex\\label{%s}'%label + text)

        exs = '\n'.join(exs)

        doc = dedent('''
        \\begin{exe}

        %s
        \\end{exe}
        ''') % exs
        
        fname = self.get_filename(
            name, 
            self.get_subdir('examples'),
            suffix='tex'
        )
        fname.write_text(doc)
        
        return df[[refcol]+textcols]