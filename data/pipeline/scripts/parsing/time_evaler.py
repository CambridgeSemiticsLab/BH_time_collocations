import json
import nbformat
import dictdiffer
from tf.fabric import Fabric
from tf.app import use
from pprint import pformat
from tools.load_parse import ParseLoader
from pathlib import Path

def mdcell(content):
    """Markdown cell for Notebook."""
    return nbformat.v4.new_markdown_cell(content)

def cdcell(code):
    """Code cell for Notebook."""
    return nbformat.v4.new_code_cell(code)

def colored_div(*spans, color=''):
    """Produce an HTML div box with color."""
    spans = '<br>'.join(str(item) for item in spans)
    div = (
        f'<div style="background:{color}; display:inline-block">'
        + spans
        + '</div>'
    )
    return div

def build_eval_notebook(todos, timeparses, paths):
    """Compile a notebook for manual review / corrections."""

    # load BHSA source_data with TF
    TF = Fabric(paths['bhsadata'], silent='deep')
    features = ( 
        'gloss function number '
        'pdp vs vt nu language '
        'rela typ label code '
    )   
    API = TF.load(features, silent='deep')
    bhsa = use('bhsa', api=API, silent='deep')
    bhsa._browse = True # ensures API.pretty outputs HTML strings
    T, F, L = bhsa.api.T, bhsa.api.F, bhsa.api.L

    # load translation source_data
    with open(paths['translations'], 'r') as infile:
        transs = json.load(infile)

    # instantiate a notebook and populate it 
    # with HTML and parsings for manual eval
    notebook = nbformat.v4.new_notebook()

    # write top of the notebook
    notebook['cells'].append(mdcell(
        f'# Time Adverbial Eval (n={len(todos)})'
    ))
    
    # format input files for the notebook
    corrs, parses = Path(paths['corrections']), Path(paths['parsed'])
    corrs, parses = corrs.resolve(), parses.resolve()

    notebook['cells'].append(cdcell(
        'from eval_tools import load_style, Tracker\n'
        'load_style()\n'
        'tracker = Tracker(\n'
        f'  "{parses}",\n'
        f'  "{corrs}",\n'
        ')\n'
        'save = tracker.save'
    ))

    for cldata in todos:

        # get needed source_data
        clause = cldata['clause']
        parsing = timeparses[clause]
        ref = str(T.sectionFromNode(clause)) 

        # build markdown cell
        cl_html = bhsa.pretty(
            clause,
            withNodes=True,
            highlights=parsing['slots'],
            hiddenTypes={'subphrase', 'clause_atom'},
        )
        esv = transs['esv'].get(ref, 'NOT AVAILABLE')
        if 'badcorr' in cldata:
            warning = colored_div(
                f'irreconcilable diffs for {clause}',
                *cldata['badcorr'],
                color='#D27988'
            ) 
        else:
            warning = ''
        markdown = mdcell(
            cl_html
            + f'ESV: {esv}<br>' 
            + warning
        )
 
        # build code cell
        indentedparse = pformat(parsing)
        code = cdcell(
            f'save({clause},\n\n'
            + f'{indentedparse}\n\n'
            + ')'
        )

        notebook['cells'].extend([markdown, code])

    # done! write the notebook to disk
    nbformat.write(notebook, paths['todo']) 

def time_evaler(paths):
    """Select and export time adverbials for eval. Apply corrections."""
    
    # load original parsings and corrections source_data
    timeparses = ParseLoader(paths['parsed']).load()
    corrections = ParseLoader(paths['corrections']).load()
    ntoget = paths['ntoget'] # number of uncorrected phrases to grab

    # sort through parsed phrases and apply corrections as found;
    # build a todo list of uncorrected parsings
    tocorrect = []
    tocorrectagain = []
    corrected = {}
    for clause, timeparse in timeparses.items():
        if clause in corrections:
            diffs = corrections[clause]
            try:
                corrected[clause] = dictdiffer.patch(diffs, timeparse)
            # underlying source_data has changed irreconcilably
            # re-issue the clause for correction
            except (KeyError, IndexError, TypeError): 
                tocorrectagain.append({'clause':clause, 'badcorr':diffs})
        else:
           tocorrect.append({'clause': clause}) 

    # export a notebook with requested number of 
    # corrections to do
    todo = (tocorrectagain + tocorrect)[:ntoget]
    build_eval_notebook(
        todo,
        timeparses,
        paths,
    )

    # export the corrected source_data
    with open(paths['evaled'], 'w') as outfile:
        json.dump(corrected, outfile, indent=2, ensure_ascii=False)
