"""
In this module, a dataset is built which
will form the basis of the quantitative
analyses.
"""

import cx_analysis.graph_nav as nav
from cx_analysis.cx import Construction
import collections
import pandas as pd
import numpy as np
from tf_tools.tokenizers import tokenize_surface
from tf_tools.formatting import book2sbl

def build_dataset(cxs, tf_api):
    """Tag features across Construction objects

    This module iterates through all Construction
    objects to make observations on each one.

    Args:
        cxs: a set of Construction objects
        tf: an instance of Text-Fabric
    """

    # set up TF shortform methods
    F, E, T, L = tf_api.F, tf_api.E, tf_api.T, tf_api.L

    # row data / observation go into here
    dataset = []

    # iterate through time adverbials and make observations
    for cx in cxs:
        
        # remove non-time adverbial phrases
        if {'not_single'} & set(cx.classification):
            continue

        # head features
        head = nav.get_headword(cx)
        head_pos = F.sp.v(head)
        head_pl = F.nu.v(head) == 'pl'
        head_du = F.nu.v(head) == 'du'
        head_sffx = F.prs.v(head) not in {'absent', 'n/a'}
        
        # phrase features
        phr_type = cx.name
        tokenized = tokenize_surface(cx.slots, tf_api, feature='lex_utf8')

        # preps
        prepositions = cx.key_roles.get('prepositions', [])
        leading_prep = next(iter(prepositions), 0)
        trailing_prep = next(iter(reversed(prepositions)), 0)
        tokenized_prep = '.'.join(F.lex_utf8.v(p) for p in prepositions)
        extended_prep = (
            not {F.pdp.v(p), F.ls.v(p)} & {'prep'}
                for p in prepositions
        )

        # quants
        quant = cx.key_roles.get('quantifier', Construction())
        quantified = 'quantified' in cx.classification
        cardinal = 'cardinal' in cx.classification
        qualitative = 'qualitative' in cx.classification
        quantifier = T.text(quant.slots) 
        quant_token = tokenize_surface(quant.slots, tf_api, feature='lex_utf8')
        
        # others
        demonstrative = cx.key_roles.get('demonstrative', 0)
        ordinal = cx.key_roles.get('ordinal', 0)
        
        # clause features
        clause = L.u(head,'clause')[0]
        cl_kind = F.kind.v(clause)
        verb = next((w for w in L.d(clause,'word') if F.pdp.v(w) == 'verb'), 0)
        tense = {'ptca':'ptcp'}.get(F.vt.v(verb), F.vt.v(verb))

        # reference features
        book,chapter,verse = T.sectionFromNode(head)
        sbl_book = book2sbl[book]
        ref = f'{sbl_book} {chapter}:{verse}'

        # set dataset defaults
        na = False #  none value
        trans = F.lex_utf8.v # transcription

        dataset.append({
            'node': L.u(head,'timephrase')[0],
            'ref': ref,
            'ph_type': phr_type,
            'text': T.text(cx.slots),
            'token': tokenized,
            'clause': T.text(clause),
            'classi': '.'.join(cx.classification),
            'time': trans(head),
            'time_etcbc': F.lex.v(head),
            'time_pos': head_pos,
            'time_pl': head_pl,
            'time_sffx': head_sffx,
            'leading_prep': trans(leading_prep) or na,
            'trailing_prep': trans(trailing_prep) or na,
            'tokenized_prep': tokenized_prep or na,
            'extended_prep': any(extended_prep) or na,
            'quantified': quantified or head_du or na,
            'quantifier': quantifier or na,
            'cardinal': quantified and cardinal, 
            'qual_quant': quant_token if (quantified and qualitative) else na,
            'demonstrative': trans(demonstrative) or na,
            'ordinal': trans(ordinal) or na,
            'cl_kind': cl_kind,
            'verb_lex': trans(verb) or na,
            'tense': tense or na,
        })

    return pd.DataFrame(dataset).set_index('node')
