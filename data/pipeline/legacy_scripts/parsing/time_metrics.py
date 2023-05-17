import re
import json
import html
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from data.pipeline.legacy_scripts.tools.html_docs import HtmlReport
from data.pipeline.legacy_scripts.tools.load_parse import ParseLoader
from pprint import pformat

def get_node_data(node, tf_api):
    """Retrieve source_data on parsed phrase for analysis."""
    F, E, L, = tf_api.F, tf_api.E, tf_api.L
    words = L.d(node, 'word')
    phatoms = L.d(node, 'phrase_atom')
    return {
        'nwords': len(words),
        'n_phatoms': len(phatoms),
    }

def build_row_data(node, tf_api, time_parsing={}, **features):
    """Build source_data that can be analyzed to assess quality of parses."""
    F, T, L = tf_api.F, tf_api.T, tf_api.L
    book, chapter, verse = T.sectionFromNode(node)

    # get tokens from error report
    err = features.get('error', '') 
    toks = re.findall('[A-Z_Ø<>/=0-9]+', err)
    tokens = (' '.join(toks)).replace('<', '(').replace('>', ')')
    if err:
        features['error'] = html.escape(err)

    # build basic source_data
    data = dict(
        node=node,
        ref=f'{book} {chapter}:{verse}',
        book=book,
        typ=F.typ.v(node),
        txt=T.text(node, fmt='text-orig-plain'),
        tokens=tokens,
        parsed=bool(time_parsing) * 1,
        time_parsing=time_parsing,
        **features
    )

    # update with time parsing source_data
    data['function'] = time_parsing.get('functions', [None])[0]
    data['slots'] = time_parsing.get('slots', [])

    # update with source_data about the phrase
    data.update(
        get_node_data(
            node, 
            tf_api
        )
    )

    return data

def build_parse_table(paths, API):
    """Produce metrics and reports on the parsings."""

    # load phrase source_data
    with open(paths['parsed'], 'r') as infile:
        parsed = json.load(infile)
    with open(paths['notparsed'], 'r') as infile:
        errors = json.load(infile)
    ph_parsings = ParseLoader(paths['ph_parses']).load()

    # build row source_data for dataframe
    rows = []
    for node, parse in parsed.items():
        node = int(node)
        row = build_row_data(
            node, 
            API, 
            time_parsing=parse
        )
        rows.append(row)
    for node, reason in errors.items():
        node = int(node)
        row = build_row_data(
            node, 
            API,
            error=str(reason)
        )
        rows.append(row)

    df = pd.DataFrame(rows)\
        .set_index('node', drop=True)
    return df

def examine_times(paths, bhsa):
    """Build a report."""

    API = bhsa.api

    # function for getting plot paths
    plotdir = Path(paths['plotsdir'])
    plotpath = lambda fname: str(plotdir.joinpath(fname))
    
    # build dataset
    df = build_parse_table(
        paths,
        API,
    )

    doc1 = HtmlReport(paths['styles'])
    doc1.heading('Time Parsing Report', 1)
    
    doc1.heading('top of table', 2)
    doc1.table(df.head())

    # Number of phrases parsed
    parsed_ct = df.parsed.value_counts()
    doc1.heading('parsed counts', 2)
    doc1.table(parsed_ct)
    parsed_pr = parsed_ct / parsed_ct.sum()
    doc1.heading('parsed perc.', 2)
    doc1.table(parsed_pr.round(2))

    # Counts of parsed semantic classes
    function_ct = df.function.value_counts()
    doc1.heading('function counts', 2)
    doc1.table(function_ct)
    function_pr = function_ct / function_ct.sum()
    doc1.heading('function perc.', 2)
    doc1.table(function_pr.round(2))

    # phrase type ct
    typ_ct = pd.pivot_table(
        df,
        index='parsed',
        columns='typ',
        aggfunc='size',
        fill_value=0,
    )
    typ_pr = typ_ct.div(typ_ct.sum(0), 1)
    doc1.heading('phrase type: parsed vs. not parsed', 2)
    doc1.table(typ_ct)
    doc1.heading('pr', 3)
    doc1.table(typ_pr)

    # plot % parsed by number of words in phrase
    nwd_ct = pd.pivot_table(
        df,
        index='nwords',
        columns='parsed',
        aggfunc='size',
        fill_value=0,
    )
    nwd_pr = nwd_ct.div(nwd_ct.sum(1), 0)    

    doc1.heading('Parsed by number of words', 2)
    doc1.heading('ct (cutoff at 20)', 3)
    doc1.table(nwd_ct[nwd_ct.index < 20].T)
    doc1.heading('pr', 3)
    doc1.table(nwd_pr[nwd_pr.index < 20].round(2).T)

    # make plot of number of words 
    fig, ax = plt.subplots(figsize=(4, 4))
    data = nwd_pr[nwd_pr.index < 20]
    x = data.index
    y = data[1]
    ax.scatter(x, y)
    ax.set_title('Number of words by % parsed (<20 words)')
    ax.set_xticks(x)
    nwd_plotpth = plotpath('nwords_by_parsed_time.svg')
    plt.savefig(nwd_plotpth, format='svg')
    doc1.img(nwd_plotpth)
   
    # counts of phrase strings that are unparsed
    err_df = df[df.parsed == 0]
    top_err_toks = err_df.tokens.value_counts().head(100)
    if top_err_toks.max() == 1:
        sortindex = top_err_toks.index.sort_values()
        top_err_toks = top_err_toks[sortindex]

    doc1.heading('most missed values', 2)
    doc1.table(top_err_toks)

    # display unparsed phrases
    doc1.table(parsed_ct)
    for string in top_err_toks.index:
        doc1.heading(string, 3)
        examples = err_df[err_df.tokens == string]
        sample_size = min([5, examples.shape[0]])
        examples = examples.sample(sample_size, random_state=42)
        for i in examples.index:
            doc1.append(
                bhsa.pretty(
                    i,
                    extraFeatures='st lex pdp function',
                    withNodes=True,
                    hiddenTypes={'subphrase', 'clause_atom'},
                )
            )
            doc1.append(err_df.loc[i]['error'])
            doc1.append('<hr>')   

    # -- BUILD PARSED INSPECTION DOCUMENT --
    doc2 = HtmlReport(paths['styles'])
    doc2.heading('Times Analysis', 1)

    # sample unique phrase temporal qualities
    for funct in function_ct.index:
        semdf = df[df.function == funct]
        size = semdf.shape[0]
        samp = semdf.sample(min(50, size), random_state=42)
        doc2.heading(funct, 3)
        for cl in samp.index:
            #function = df.loc[cl]['function']
            parsing = pformat(df.loc[cl]['time_parsing'])
            doc2.append(
                bhsa.plain(
                    cl, 
                ).replace('div class="rtl"', 'div')
            )   
            doc2.append(f'{cl}<br>')
            doc2.append(html.escape(f'{parsing}'))
            doc2.append('<br><br>')
        doc2.append('<hr>')

    # finish and export
    doc1.export()
    with open(paths['metrics'], 'w') as outfile:
        outfile.write(doc1.html)

    # finish and export
    doc2.export()
    with open(paths['catmetrics'], 'w') as outfile:
        outfile.write(doc2.html)
