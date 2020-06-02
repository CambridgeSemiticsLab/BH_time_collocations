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

class BoolMap:
    """Set and map boolean values"""
    def __init__(self, true=True, false=False):
        self.true = true
        self.false = false

    def map(expression): 
        """Convert boolean to mapped values.""" 
        if expression:
            return self.true
        else:
            return self.false 

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

    # set dataset defaults
    true = True
    false = False
    null = np.nan

    # function to map expressions
    def boomap(expression):
        if expression:
            return true
        else:
            return false

    trans = F.lex_utf8.v # transcription

    # iterate through time adverbials and make observations
    for cx in cxs:
        
        # avoid non-time adverbial phrases
        if {'not_single'} & set(cx.classification):
            continue

        # cx
        cxclass = cx.classification

        # head features
        head = nav.get_headword(cx)
        head_du = F.nu.v(head) == 'du'
        
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
        quantified = 'quantified' in cxclass
        qualitative = 'qualitative' in cxclass
        quant_str = T.text(quant.slots) 
        quant_token = tokenize_surface(quant.slots, tf_api, feature='lex_utf8')
        
        # others
        demon_cx = cx.key_roles.get('demonstrative', 0)
        ordinal_cx = cx.key_roles.get('ordinal', 0)
        
        # clause features
        clause = L.u(head,'clause')[0]
        sentence = L.u(head,'sentence')[0]
        verb = next((w for w in L.d(clause,'word') if F.pdp.v(w) == 'verb'), 0)
        tense = {'ptca':'ptcp'}.get(F.vt.v(verb), F.vt.v(verb))

        # reference features
        book,chapter,verse = T.sectionFromNode(head)
        sbl_book = book2sbl[book]
        ref = f'{sbl_book} {chapter}:{verse}'

        # phrase
        phrase = L.u(head,'phrase')[0]
        function = F.function.v(phrase)

        # mappings of features
        # map demonstratives to near/far labels
        demon_map = {
            'Z>T': 'near',
            'HJ>': 'far',
            'HMH': 'far',
            '>LH': 'near',
            'HM': 'far',
            'HW>': 'far',
            'ZH': 'near'
        }        
        data = {
            'node': phrase,
            'function': function,
            'ref': ref,
            'book': book,
            'ph_type': phr_type,
            'head': trans(head),
            'text': T.text(cx.slots),
            'token': tokenized,
            'clause': T.text(clause),
            'sentence': T.text(sentence),
            'classi': '.'.join(cxclass),
            'head_etcbc': F.lex.v(head),
            'head_pos': F.sp.v(head),
            'plural': boomap(F.nu.v(head) == 'pl'),
            'suffix': boomap(F.prs.v(head) not in {'absent', 'n/a'}), 
            'preposition': boomap('prep' in cxclass),
            'leading_prep': trans(leading_prep) or null,
            'trailing_prep': trans(trailing_prep) or null,
            'tokenized_prep': tokenized_prep or null,
            'extended_prep': boomap(any(extended_prep)),
            'ø': boomap('bare' in cxclass),
            'øanchor': boomap('øanchor' in cxclass),
            'genitive': boomap('genitive' in cxclass),
            'definite': boomap('definite' in cxclass),
            'quantified': boomap(quantified or head_du), 
            'quant_str': quant_str or null,
            'cardinal': boomap(quantified and 'cardinal' in cxclass), 
            'qualitative': boomap(qualitative),
            'qual_str': quant_token if (quantified and qualitative) else null,
            'demonstrative': boomap('demonstrative' in cxclass),
            'demon_str': trans(demon_cx) or null,
            'demon_dist': demon_map.get(F.lex.v(demon_cx)) or null,
            'ordinal': boomap('ordinal' in cxclass),
            'ord_str': trans(ordinal_cx) or null,
            'cl_kind': F.kind.v(clause),
            'verb': boomap(bool(verb)),
            'tense': tense or null,
            'verb_lex': trans(verb) or null,
            'book_sbl': sbl_book,
            'lang': F.language.v(head),
        }
    
        dataset.append(data)

    return pd.DataFrame(dataset).set_index(['node'])
