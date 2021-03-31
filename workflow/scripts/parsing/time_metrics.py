import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tools.html_docs import HtmlReport
import tools.nav_tree as nt
from tools.load_parse import ParseLoader

def build_row_data(node, tf_api, ph_parse, time_parsing=None, **features):
    """Build data that can be analyzed to assess quality of parses."""
    F, T, L = tf_api.F, tf_api.T, tf_api.L
    book, chapter, verse = T.sectionFromNode(node)
    phrase = L.u(node, 'phrase')[0]
    words = L.d(node, 'word')
    if len(ph_parse) > 1:
        head = nt.get_head(ph_parse)
    else:
        head = ph_parse[0]
    return dict(
        node=node,
        ref=f'{book} {chapter}:{verse}',
        book=book,
        parse=time_parsing,
        parsed=bool(time_parsing) * 1,
        typ=F.typ.v(node),
        nwords=len(words),
        txt=T.text(node, fmt='text-orig-plain'),
        head=T.text(head),
        head_lex=F.lex.v(head),
        **features
    )

def build_parse_table(parsed_path, notparsed_path, ph_parsings, API):
    """Produce metrics and reports on the parsings."""

    # load the parsed/not parsed files
    with open(parsed_path, 'r') as infile:
        parsed = json.load(infile)
    with open(notparsed_path, 'r') as infile:
        errors = json.load(infile)

    # build row data for dataframe
    rows = []
    for node, parse in parsed.items():
        node = int(node)
        row = build_row_data(
            node, 
            API, 
            ph_parsings[node],
            time_parsing=parse
        )
        rows.append(row)
    for node, reason in errors.items():
        node = int(node)
        row = build_row_data(
            node, 
            API,
            ph_parsings[node],
            error=str(reason)
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    return df

def examine_times(parsed_path, notparsed_path, 
    datalocs, stylepaths, metricspath, bhsa=None):
    """Build a report."""

    API = bhsa.api

    # load custom pos tags
    ph_parses = ParseLoader(datalocs['ph_parsings']).load()

    # function for getting plot paths
    plotdir = Path(datalocs['plotsdir'])
    plotpath = lambda fname: str(plotdir.joinpath(fname))
    
    # build dataset
    df = build_parse_table(
        parsed_path, 
        notparsed_path,
        ph_parses,
        API,
    )
    
    doc = HtmlReport(stylepaths)
    doc.heading('Time Parsing Report', 1)
    
    doc.heading('top of table', 2)
    doc.table(df.head())

    # Number of phrases parsed
    parsed_ct = df.parsed.value_counts()
    doc.heading('parsed counts', 2)
    doc.table(parsed_ct)
    parsed_pr = parsed_ct / parsed_ct.sum()
    doc.heading('parsed perc.', 2)
    doc.table(parsed_pr.round(2))

    # Counts of parsed semantic classes
    semclass_ct = df.parse.value_counts()
    doc.heading('semclass counts', 2)
    doc.table(semclass_ct)
    semclass_pr = semclass_ct / semclass_ct.sum()
    doc.heading('semclass perc.', 2)
    doc.table(semclass_pr.round(2))

    # Number of phrases parsed/not parsed by type
    funct_ct = pd.pivot_table(
        df,
        index='parsed',
        columns='typ',
        aggfunc='size',
        fill_value=0,
    )
    funct_pr = funct_ct.div(funct_ct.sum(0), 1)
    doc.heading('phrase type: parsed vs. not parsed', 2)
    doc.table(funct_ct)
    doc.heading('pr', 3)
    doc.table(funct_pr)

   # plot % parsed by number of words in phrase
    nwd_ct = pd.pivot_table(
        df,
        index='nwords',
        columns='parsed',
        aggfunc='size',
        fill_value=0,
    )
    nwd_pr = nwd_ct.div(nwd_ct.sum(1), 0)    

    doc.heading('Parsed by number of words', 2)
    doc.heading('ct (cutoff at 20)', 3)
    doc.table(nwd_ct[nwd_ct.index < 20].T)
    doc.heading('pr', 3)
    doc.table(nwd_pr[nwd_pr.index < 20].round(2).T)

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
    doc.img(nwd_plotpth)

    # look at DP for books
#    book_ct = pd.pivot_table(
#        df,
#        index='book',
#        columns=['function', 'parsed'],
#    )
    
    # counts of phrase strings that are unparsed
    err_df = df[df.parsed == 0]
    top_err_str = err_df.error.value_counts().head(50)
    doc.heading('most missed values', 2)
    doc.table(top_err_str)

    # display unparsed phrases
    for string in top_err_str.index:
        doc.heading(string, 3)
        examples = err_df[err_df.error == string]
        sample_size = min([5, examples.shape[0]])
        examples = examples.sample(sample_size, random_state=42)
        for i in examples.index:
            doc.append(
                bhsa.pretty(
                    err_df.loc[i]['node'], 
                    extraFeatures='typ st',
                    withNodes=True
                )
            )
            doc.append(err_df.loc[i]['error'])
        doc.append('<hr>')

    # finish and export
    doc.export()
    with open(metricspath, 'w') as outfile:
        outfile.write(doc.html)
