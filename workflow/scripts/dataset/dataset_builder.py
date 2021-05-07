import pandas as pd
from tf.fabric import Fabric
from tools.load_parse import ParseLoader
from bidi.algorithm import get_display

# custom modules
from .tokenizers import tokenize_lexemes
from .nav_time import get_times
from .verb_form import get_verbform
from .modis_getter import get_modis
from .synvar_carc import in_dep_calc as clause_relator
from .modify_domain import permissive_q
from .book_formats import get_book_maps, etcbc2sbl, etcbc2abbr

# logic for parsing the features
function_simp = {
    'cal_simul': 'simultaneous',
    'multi_simul': 'simultaneous',
    'multi_habitual': 'habitual',
}

quality_map = {
    'simultaneous': 'location',
    'atelic_ext': 'duration',
    'anterior_dur': 'duration',
    'posterior': 'sequence',
    'posterior_dur': 'duration',
    'habitual': 'iteration',
    'begin_to_end': 'duration',
#    'purposive_ext': 'duration',
    'simultaneous + atelic_ext': 'duration',
    'regular_recurrence': 'iteration',
    'multi_simuls': 'iteration',
    'anterior': 'sequence',
#    'telic_ext': 'duration',
#    'dist_posterior': 'duration',
#    'dist_future': 'duration',
    'simul_to_end': 'duration',
    'begin_to_end_habitual': 'iteration',
}

def get_clause_data(clause, API):
    """Build data on a clause."""
    F, E, T, L = API.F, API.E, API.T, API.L
    cl_rela = clause_relator(clause, API)
    verbs = [
        w for w in L.d(clause,'word')
            if F.pdp.v(w) == 'verb'
    ]
    verb = verbs[0] if verbs else None
    verb_lex = F.lex.v(verb) if verb else None
    verb_txt = F.lex_utf8.v(verb)
    if verb:
        verb_form = get_verbform(verb, API)
    else:
        verb_form = None
    data = {
        'cl_rela': cl_rela,
        'verb': verb,
        'verbform': verb_form,
        'verb_etcbc': verb_lex,
        'verb_txt': verb_txt,
    }
    return data

def build_dataset(paths):
    """Construct tabular dataset of time adverbials."""

    # load needed TF BHSA data
    TF = Fabric(locations=paths['bhsadata'], silent='deep')
    API = TF.load(
        'kind lex vt pdp ls '
        'lex_utf8 nu code rela '
        'prs prs_gn prs_nu prs_ps '
        'genre mother txt uvf '
    )
    F, E, T, L = API.F, API.E, API.T, API.L

    # get map of books to various subcollections
    bookmap = get_book_maps(API)

    # load data parsed in this project
    time_data = ParseLoader(paths['timedata']).load()
    phrase_data = ParseLoader(paths['phrasedata']).load()

    # build the rows  
    rows = []
    for clause, data in time_data.items():
        slots = sorted(data['slots'])
        phrases = data['phrase_nodes']
        ph_parses = [phrase_data[ph]['parse'] for ph in phrases]
        text = T.text(slots)
        versen = L.u(slots[0], 'verse')[0]
        genre = F.genre.v(versen)
        domain = permissive_q(clause, API)
        gendom = f'{genre}.{domain}'
        cl_text = T.text(clause)
        sentn = L.u(clause, 'sentence')[0]
        sent_text = T.text(sentn)
        bk, ch, vs = T.sectionFromNode(versen)
        book_super = bookmap['super'].get(bk, bk)
        canon_part = bookmap['tripart'][bk]
        period = bookmap['period'].get(bk, '')
        verse = f'{bk} {ch}:{vs}'
        name = data['functions'][0]
        function = function_simp.get(name, name)
        quality = quality_map.get(function, None)
        times = sorted(get_times(data))
        times_etcbc = '|'.join(F.lex.v(t) for t in times)
        times_utf8 = '|'.join(F.lex_utf8.v(t) for t in times)
        times_utf8d = (
            get_display(times_utf8)
                .replace('\u05c1', '') # rm shin dots b/c don't display well
                .replace('\u05c2', '')
        )
        lex_token = tokenize_lexemes(
            ph_parses, 
            API,
            heads=times,
        )
        if len(times) == 1:
            is_advb = 1*(F.pdp.v(times[0]) == 'advb')
        else:
            is_advb = 0 
  
        data = {
            'node': clause,
            'verse': verse,
            'book': bk,
            'booksuper': book_super,
            'canon_part': canon_part,
            'period': period,
            'genre': genre,
            'domain': domain,
            'gendom': gendom,
            'function': function,
            'quality': quality,
            'name': name,
            'text': text,
            'times_etcbc': times_etcbc,
            'times_utf8': times_utf8,
            'times_utf8d': times_utf8d,
            'n_times': len(times),
            'lex_token': lex_token,
            'is_advb': is_advb,
            'cl_kind': F.kind.v(clause),
            'clause': cl_text,
            'sentence': sent_text,
        }

        # add clause-based data
        data.update(
            get_clause_data(clause, API)
        )
        
        # add modifier data
        # this data will typically be used to analyze single-phrased TAs
        # so we just take the first ph_parse
        modifiers = get_modis(ph_parses[0], API)
        data.update(
            modifiers
        )
    
        # mark unmodified words as adverbs
        if not modifiers or (len(modifiers) == 1 and 'PP' in modifiers):
            data['unmodified'] = 1

        # finish
        rows.append(data)

    # turn into data frame and export as CSV
    df = pd.DataFrame(rows) 
    print('exporting:', df.shape)
    df.to_csv(paths['dataset'], index=False)
