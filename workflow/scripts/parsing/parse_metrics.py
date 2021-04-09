import re
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tf.app import use
from tf.fabric import Fabric
from pathlib import Path
from tools import nav_tree as nt

# NB custom module
from tools.html_docs import HtmlReport

def build_row_data(node, tf_api, slot2pos, parsing=[], **features):
    """Build data that can be analyzed to assess quality of parses."""
    F, T, L = tf_api.F, tf_api.T, tf_api.L
    book, chapter, verse = T.sectionFromNode(node)
    phrase = L.u(node, 'phrase')[0]
    words = L.d(node, 'word')

    # get tokens from error report
    err = features.get('error', '')
    toks = re.findall('[A-Z_]+', err)
    tokens = ' '.join(toks)

    return dict(
        node=node,
        ref=f'{book} {chapter}:{verse}',
        book=book,
        parse=str(parsing),
        parsed=bool(parsing) * 1,
        typ=F.typ.v(node),
        function=F.function.v(phrase),
        nwords=len(words),
        txt=T.text(node, fmt='text-orig-plain'),
        tokens=tokens,
        **features
    )

def build_parse_tables(paths, API):
    """Produce metrics and reports on the parsings."""
    
    # load phrase and POS data
    with open(paths['parsed_atoms'], 'r') as infile:
        parsed = json.load(infile)
    with open(paths['unparsed_atoms'], 'r') as infile:
        errors = json.load(infile)
    with open(paths['slot2pos'], 'r') as infile:
        slot2pos = json.load(infile)

    rows = []
    sp_rows = []
    sp_id = 1 # counter to make subphrase ids

    # build data for successful parses
    for node, parse in parsed.items():
        node = int(node)
        row = build_row_data(
            node, 
            API, 
            slot2pos,
            parsing=parse
        )
        rows.append(row)

        # add in subphrase data
        if type(parse) == int or len(parse) < 3:
            continue
        subphrases = list(nt.traverse_tree(parse))
        for sp in subphrases:
            src, tgt, kind = sp
            slots = list(nt.get_slots(sp))
            nslots = len(slots)
            sp_rows.append({
                'id': sp_id,
                'phrase': node,
                'nslots': len(slots),
                'kind': kind,
                'parse': sp,
            })
            sp_id += 1

    # build data for unsuccesful parses
    for node, reason in errors.items():
        node = int(node)
        row = build_row_data(
            node, 
            API,
            slot2pos,
            error=str(reason)
        )
        rows.append(row)

    df = pd.DataFrame(rows)\
        .set_index('node', drop=True)
    sp_df = pd.DataFrame(sp_rows)\
        .set_index('id', drop=True)

    return df, sp_df

def examine_parsings(paths, bhsa):
    """Build a HTML report about phrase parsings."""
    
    API = bhsa.api

    # function for getting plot paths
    plotdir = Path(paths['plotsdir'])
    plotpath = lambda fname: str(plotdir.joinpath(fname))
    
    # build datasets for doing analysis on 
    df, sp_df = build_parse_tables(paths, API)

    # optionally restrict to Time Phrases
#    df = df[df.function == 'Time']
#    sp_df = sp_df[sp_df.phrase.isin(df.index)]

    # -- BUILD UN-PARSED INSPECTION DOCUMENT --

    doc1 = HtmlReport(paths['styles'])
    doc1.heading('Parsing Report', 1)

    doc1.heading('top of table', 2)
    doc1.table(df.head())

    # Number of phrases parsed
    parsed_ct = df.parsed.value_counts()
    doc1.heading('parsed counts', 2)
    doc1.table(parsed_ct)
    parsed_pr = parsed_ct / parsed_ct.sum()
    doc1.heading('parsed perc.', 2)
    doc1.table(parsed_pr.round(2))

    # Number of phrases parsed/not parsed by function
    funct_ct = pd.pivot_table(
        df,
        index='parsed',
        columns='function',
        aggfunc='size',
        fill_value=0,
    )
    funct_pr = funct_ct.div(funct_ct.sum(0), 1)
    doc1.heading('function parsed vs. not parsed', 2)
    doc1.table(funct_ct)
    doc1.heading('pr', 3)
    doc1.table(funct_pr)

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
    nwd_plotpth = plotpath('nwords_by_parsed.svg')
    plt.savefig(nwd_plotpth, format='svg')
    doc1.img(nwd_plotpth)

    # look at DP for books
    book_ct = pd.pivot_table(
        df,
        index='book',
        columns=['function', 'parsed'],
    )

    # counts of phrase strings that are unparsed
    err_df = df[df.parsed == 0]
    top_err_toks = err_df.tokens.value_counts().head(50)
    doc1.heading('most missed values', 2)
    doc1.table(top_err_toks)

    # display unparsed phrases
    for string in top_err_toks.index:
        doc1.heading(string, 3)
        examples = err_df[err_df.tokens == string]
        sample_size = min([5, examples.shape[0]])
        examples = examples.sample(sample_size, random_state=42)
        for i in examples.index:
            doc1.append(
                bhsa.pretty(
                    i,
                    extraFeatures='typ st',
                    withNodes=True
                )
            )
            doc1.append(err_df.loc[i]['error'])
            doc1.append('<hr>')

    # -- BUILD PARSED INSPECTION DOCUMENT --
    doc2 = HtmlReport(paths['styles'])
    doc2.heading('Parsed Phrase Analysis', 1)

    doc2.heading('Phrase type counts', 2)
    kind_ct = sp_df.kind.value_counts()
    doc2.table(kind_ct)

    # get random examples for each kind
    for kind in kind_ct.index:
        kdf = sp_df[sp_df.kind == kind]
        ksize = kdf.shape[0]
        samp = kdf.sample(min(50, ksize), random_state=42)
        doc2.heading(kind, 3)
        for i in samp.index:
            ph = kdf.loc[i]['phrase']
            sp_parse = eval(str(kdf.loc[i]['parse']))
            sp_slots = nt.get_slots(sp_parse)
            ph_parse = eval(str(df.loc[ph]['parse']))
            ph_show_parse = nt.show_relas(
                ph_parse, 
                API.T.text,
                '<br>'
            )
            doc2.append(
                bhsa.plain(
                    ph,
                    highlights=sp_slots,
                    withNodes=True,
                ).replace('div class="rtl"', 'div')
            )
            doc2.append(f'{ph}<br>')
            doc2.append(ph_show_parse) 
            doc2.append('<br><br>')
        doc2.append('<hr>')

    # finish and export documents
    doc1.export()
    with open(paths['parsedmets'], 'w') as outfile:
        outfile.write(doc1.html)

    doc2.export()
    with open(paths['notparsedmets'], 'w') as outfile:
        outfile.write(doc2.html)
