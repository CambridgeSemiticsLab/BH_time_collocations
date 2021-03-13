import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tf.app import use
from tf.fabric import Fabric
from pathlib import Path

# NB custom module
from tools.html_docs import HtmlReport

def get_pos_chain(node, tf, slot2pos):
    """Build pdp string chain for a phrase."""
    pos_chain = []
    for w in tf.L.d(node, 'word'):
        pos = slot2pos[str(w)]
        pos_chain.append(pos)
        if tf.F.st.v(w) == 'c' and pos != 'PREP':
            pos_chain.append('C')
        if tf.F.prs.v(w) not in {'absent', 'n/a'}:
            pos_chain.append('SFX')
    return ' '.join(pos_chain) 

def build_row_data(node, tf_api, slot2pos, parsing=[], **features):
    """Build data that can be analyzed to assess quality of parses."""
    F, T, L = tf_api.F, tf_api.T, tf_api.L
    book, chapter, verse = T.sectionFromNode(node)
    pos_chain = get_pos_chain(node, tf_api, slot2pos)
    phrase = L.u(node, 'phrase')[0]
    words = L.d(node, 'word')
    return dict(
        node=node,
        ref=f'{book} {chapter}:{verse}',
        parse=str(parsing),
        parsed=bool(parsing) * 1,
        typ=F.typ.v(node),
        function=F.function.v(phrase),
        nwords=len(words),
        txt=T.text(node, fmt='text-orig-plain'),
        pos_chain=pos_chain,
        **features
    )

def build_parse_table(parsed_path, notparsed_path, datalocs, API, slot2pos):
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
            slot2pos,
            parsing=parse
        )
        rows.append(row)
    for node, reason in errors.items():
        node = int(node)
        row = build_row_data(
            node, 
            API,
            slot2pos,
            error=str(reason)
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    return df

def examine_parsings(parsed_path, notparsed_path, 
    datalocs, stylepaths, metricspath, bhsa=None):
    """Build a report."""

    # initialize TF if necessary
    if bhsa is None:
        TF = Fabric(locations=datalocs['bhsadata'])
        features = (
            'rela code gloss function number '
            'pdp vs vt typ language label st '
            'prs'
        )
        API = TF.load(features, silent='deep')
        bhsa = use('bhsa', api=API, silent='deep')
        bhsa._browse = True
        API = bhsa.api
    else:
        API = bhsa.api

    # load custom pos tags
    with open(datalocs['slot2pos'], 'r') as infile:
        slot2pos = json.load(infile)

    # function for getting plot paths
    plotdir = Path(datalocs['plotsdir'])
    plotpath = lambda fname: str(plotdir.joinpath(fname))
    
    # build dataset
    df = build_parse_table(
        parsed_path, 
        notparsed_path,
        datalocs,
        API,
        slot2pos
    )
    
    rep = HtmlReport(stylepaths)
    rep.heading('Parsing Report', 1)
    
    rep.heading('top of table', 2)
    rep.table(df.head())

    # Number of phrases parsed
    parsed_ct = df.parsed.value_counts()
    rep.heading('parsed counts', 2)
    rep.table(parsed_ct)
    parsed_pr = parsed_ct / parsed_ct.sum()
    rep.heading('parsed perc.', 2)
    rep.table(parsed_pr.round(2))

    # Number of phrases parsed/not parsed by function
    funct_ct = pd.pivot_table(
        df,
        index='parsed',
        columns='function',
        aggfunc='size',
        fill_value=0,
    )
    funct_pr = funct_ct.div(funct_ct.sum(0), 1)
    rep.heading('function parsed vs. not parsed', 2)
    rep.table(funct_ct)
    rep.heading('pr', 3)
    rep.table(funct_pr)
   
    # plot % parsed by number of words in phrase
    nwd_ct = pd.pivot_table(
        df,
        index='nwords',
        columns='parsed',
        aggfunc='size',
        fill_value=0,
    )
    nwd_pr = nwd_ct.div(nwd_ct.sum(1), 0)    
    fig, ax = plt.subplots(figsize=(4, 4))
    data = nwd_pr[nwd_pr.index < 20]
    x = data.index
    y = data[1]
    ax.scatter(x, y)
    ax.set_title('Number of words by % parsed (<20 words)')
    ax.set_xticks(x)
    nwd_plotpth = plotpath('nwords_by_parsed.svg')
    plt.savefig(nwd_plotpth, format='svg')
    rep.img(nwd_plotpth)

    # counts of phrase strings that are unparsed
    err_df = df[df.parsed == 0]
    top_err_str = err_df.pos_chain.value_counts().head(25)
    rep.heading('most missed values', 2)
    rep.table(top_err_str)

    # display unparsed phrases
    for string in top_err_str.index:
        rep.heading(string, 3)
        examples = err_df[err_df.pos_chain == string]
        examples = examples.sample(5, random_state=42)
        for i in examples.index:
            rep.append(
                bhsa.pretty(
                    err_df.loc[i]['node'], 
                    extraFeatures='typ st',
                    withNodes=True
                )
            )
            rep.append(err_df.loc[i]['error'])
        rep.append('<hr>')

    # finish and export
    rep.export()
    with open(metricspath, 'w') as outfile:
        outfile.write(rep.html)
