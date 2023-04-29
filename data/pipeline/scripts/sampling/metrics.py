import json
import pandas as pd
import numpy as np
from tf.fabric import Fabric

# NB custom module stored at ../
from tools.html_docs import HtmlReport  

def get_rows(node, tf_api, **features):
    """Get data on a given phrase node."""
    F, T, L = tf_api.F, tf_api.T, tf_api.L
    book, chapter, verse = T.sectionFromNode(node)
    ref = f'{book} {chapter}:{verse}'
    phrase = L.u(node, 'phrase')[0]
    data = dict( 
        node=node,
        ref=ref,
        book=book,
        typ=F.typ.v(node),
        plain_txt=T.text(node, fmt='text-orig-plain'),
        nwords=len(L.d(node, 'word')),
        nsubphrases=len(L.d(node, 'subphrase')),
        parent_atoms=len(L.d(phrase, 'phrase_atom')),
    )
    data.update(**features)
    return data

def build_dataset(samp_path, nosamp_path, tf_locs):
    """Build DataFrame for analysis"""
    # initialize TF with BHSA data
    TF = Fabric(locations=tf_locs)
    API = TF.load('typ')
    F, T, L = API.F, API.T, API.L 

    # load the samples
    with open(samp_path, 'r') as infile:
        samps = json.load(infile)
    with open(nosamp_path, 'r') as infile:
        nosamps = json.load(infile)

    # build table for analysis
    rows = []
    for node in samps:
        rows.append(
            get_rows(
                int(node), 
                API,
                kept=True
            )
        )
    for reason, nodes in nosamps.items():
        for node in nodes:
            rows.append(
                get_rows(
                    node, 
                    API, 
                    reason=reason,
                    kept=False
                )
            )

    df = pd.DataFrame(rows)
    return df

def measure_phrases(samp_path, nosamp_path, 
                    tf_locs, outpath, outdir,
                    stylesheets=[]):
    """Produce metrics for sampled / non-sampled phrases."""

    df = build_dataset(
        samp_path,
        nosamp_path,
        tf_locs
    )

    # count breakdown between keep / reject overall
    keep_ct = df.kept.value_counts()
    keep_pr = keep_ct / keep_ct.sum()

    # count breakdown by phrase type
    ptype_ct = pd.pivot_table(
        df*1,
        index='kept',
        columns='typ',
        aggfunc='size',
        fill_value=0,
    )
    ptype_pr = ptype_ct.div(ptype_ct.sum(0), 1)
    ptype_dpr = ptype_pr - ptype_pr.iloc[::-1].values
    
    # look at rejection reasons
    rej_df = df[df.kept == False]
    reason_ct = rej_df.reason.value_counts()

    # build HTML
    rep = HtmlReport(stylesheets=stylesheets)
    rep.heading('Sample Stats')
    rep.heading('keep ct', 2)
    rep.table(keep_ct)
    rep.heading('keep pr', 2)
    rep.table(keep_pr)
    rep.heading('keep ct by type', 2)
    rep.table(ptype_ct.round(2))    
    rep.heading('keep pr by type', 2)
    rep.table(ptype_pr.round(3))
    rep.heading('keep dpr by type', 2)
    rep.table(ptype_dpr.round(3))
    rep.heading('rejection reasons', 2)
    rep.table(reason_ct)
    rep.export()
    with open(outpath, 'w') as outfile:
        outfile.write(rep.html)
