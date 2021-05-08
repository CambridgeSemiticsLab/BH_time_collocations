import pandas as pd
from tf.fabric import Fabric
from bidi.algorithm import get_display

# custom modules
import tools.nav_tree as nt
from tools.load_parse import ParseLoader
from .tokenizers import tokenize_lexemes
from .nav_time import get_times
from .verb_form import get_verbform
from .modis_getter import get_modis
from .synvar_carc import in_dep_calc as clause_relator
from .modify_domain import permissive_q
from .book_formats import get_book_maps, etcbc2sbl, etcbc2abbr

def remove_shindots(string):
    """Remove dots from ש"""
    return(
        string
            .replace('\u05c1', '') 
            .replace('\u05c2', '')
    )

def get_ref_data(node, API):
    """Retrieve refrence data."""

    # TF methods
    F, E, T, L = API.F, API.E, API.T, API.L

    # get map of books to various subcollections
    bookmap = get_book_maps(API)
    
    # calc various ref data
    bk, ch, vs = T.sectionFromNode(node)
    book_super = bookmap['super'].get(bk, bk)
    canon_part = bookmap['tripart'][bk]
    period = bookmap['period'].get(bk, '')
    verse = f'{bk} {ch}:{vs}'

    # package and ship
    ref_data = {
        'verse': verse,
        'book': bk,
        'booksuper': book_super,
        'canon_part': canon_part,
        'period': period,
    }
    return ref_data 

def get_clause_data(clause, API):
    """Build data on a clause."""

    # TF methods
    F, E, T, L = API.F, API.E, API.T, API.L

    # collect data of clause, but also
    # of embedding linguistic unit like sentence
    firstw = L.d(clause, 'word')[0]
    versen = L.u(firstw, 'verse')[0]
    genre = F.genre.v(versen)
    domain = permissive_q(clause, API)
    gendom = f'{genre}.{domain}'
    cl_text = T.text(clause)
    sentn = L.u(clause, 'sentence')[0]
    sent_text = T.text(sentn)
    cl_rela = clause_relator(clause, API)

    # build data for clause verbs
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

    # package and ship it!
    data = {
        'verb': verb,
        'verbform': verb_form,
        'verb_etcbc': verb_lex,
        'verb_txt': verb_txt,
        'clause': cl_text,
        'sentence': sent_text,
        'cl_rela': cl_rela,
        'domain': domain,
        'genre': genre,
        'gendom': gendom,
        'cl_kind': F.kind.v(clause),
    }
    return data

def get_word_formats(words, API, 
                     prefix='', joiner='|', 
                     latexrow='words_latex'):
    """Retrieve word formats for a list of words."""

    # BHSA methods
    F = API.F

    # format data on heads
    words_etcbc = f'{joiner}'.join(F.lex.v(t) for t in words)
    words_utf8 = f'{joiner}'.join(F.lex_utf8.v(t) for t in words)
    words_utf8d = (
        get_display(remove_shindots(words_utf8))
    )
    # formatting for inclusion in latex docs (esp. tables)
    words_latex = remove_shindots('\texthebrew{%s}' % words_utf8)

    # package and ship
    formats = {
        f'{prefix}etcbc': words_etcbc,
        f'{prefix}utf8': words_utf8,
        f'{prefix}utf8d': words_utf8d,
        f'{latexrow}': words_latex,
    }
    return formats

def time_dataset(paths, API):
    """Construct tabular dataset of time adverbials."""

    print('Building time dataset...')

    # text-fabric methods
    F, E, T, L = API.F, API.E, API.T, API.L

    # logic for parsing time features
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
        'simultaneous + atelic_ext': 'duration',
        'regular_recurrence': 'iteration',
        'multi_simuls': 'iteration',
        'anterior': 'sequence',
        'simul_to_end': 'duration',
        'begin_to_end_habitual': 'iteration',
    #    'purposive_ext': 'duration',
    #    'telic_ext': 'duration',
    #    'dist_posterior': 'duration',
    #    'dist_future': 'duration',
     
    }

    # load data parsed in this project
    time_data = ParseLoader(paths['timedata']).load()
    phrase_data = ParseLoader(paths['phrasedata']).load()

    # build the rows  
    rows = []
    modi_keys = set()
    for clause, data in time_data.items():

        # build up the row data here
        rowdata = {
            'node': clause,
        }

        # add reference data
        rowdata.update(
            get_ref_data(clause, API)
        )

        # add formatted, | separated lexeme strings for timewords
        times = sorted(get_times(data))
        rowdata.update(
            get_word_formats(
                times, 
                API, 
                prefix='times_',
                latexrow='TA Heads',
            )
        ) 

        # add data related to the phrase
        slots = sorted(data['slots'])
        phrases = data['phrase_nodes']
        ph_parses = [phrase_data[ph]['parse'] for ph in phrases]
        text = T.text(slots)
        name = data['functions'][0]
        function = function_simp.get(name, name)
        quality = quality_map.get(function, None)
        lex_token = tokenize_lexemes(
            ph_parses, 
            API,
            heads=times,
        )
        if len(times) == 1:
            is_advb = 1*(F.pdp.v(times[0]) == 'advb')
        else:
            is_advb = 0 

        rowdata.update({
            'function': function,
            'quality': quality,
            'name': name,
            'text': text,
            'n_times': len(times),
            'lex_token': lex_token,
            'is_advb': is_advb,
        })

        # add clause-based data
        rowdata.update(
            get_clause_data(clause, API)
        )
        
        # add modifier data
        # this data will typically be used to analyze single-phrased TAs
        # so we just take the first ph_parse
        modifiers = get_modis(ph_parses[0], API)
        modi_keys |= set(modifiers)
        rowdata.update(
            modifiers
        )
    
        # mark unmodified words as adverbs
        if not modifiers or (len(modifiers) == 1 and 'PP' in modifiers):
            rowdata['unmodified'] = 1

        # finish
        rows.append(rowdata)

    # turn into data frame and export as CSV
    df = pd.DataFrame(rows) 
    modi_keys = list(modi_keys)
    df.loc[:, modi_keys] = df.loc[:, modi_keys].fillna(0)
    print('\texporting time df:', df.shape)
    df.to_csv(paths['timedataset'], index=False)

def phrase_dataset(paths, API):
    """Construct tabular dataset of time adverbials."""

    print('Building phrase dataset...')

    # text-fabric methods
    F, E, T, L = API.F, API.E, API.T, API.L

    # load data parsed in this project
    phrase_data = ParseLoader(paths['phrasedata']).load()
    functions = ParseLoader(paths['functions']).load()

    # build the rows  
    rows = []
    modi_keys = set()
    for phrase, ph_data in phrase_data.items():

        # build up the row data here
        rowdata = {
            'node': phrase,
        }

        # add reference data
        rowdata.update(
            get_ref_data(phrase, API)
        )

        # load / configure parsing data
        parse = ph_data['parse']
        if len(parse) == 1:
            parse = [None, parse[0], None]

        # add formatted, | separated lexeme strings for timewords
        subphrases = list(nt.unfold_paras(parse))
        heads = [
            nt.get_head(sp) for sp in subphrases
        ]
        rowdata.update(
            get_word_formats(
                heads, 
                API, 
                prefix='heads_',
                latexrow='Phrase Heads',
            )
        ) 

        # add data related to the phrase
        slots = sorted(
            s for sp in subphrases
                for s in nt.get_slots(sp)
        )
        phrases = sorted(ph_data['has_phrases'])
        head_ph = L.u(heads[0], 'phrase')[0]
        text = T.text(slots)
        function = functions[head_ph]
        rowdata.update({
            'function': function,
            'text': text,
            'n_heads': len(heads),
        })

        # add clause-based data
        clause = L.u(head_ph, 'clause')[0]
        rowdata['clause_node'] = clause
        rowdata.update(
            get_clause_data(clause, API)
        )
        
        # add modifier data
        modifiers = get_modis(parse, API)
        modi_keys |= set(modifiers)
        rowdata.update(
            modifiers
        )
    
        # mark unmodified phrase heads as such
        if not modifiers or (len(modifiers) == 1 and 'PP' in modifiers):
            rowdata['unmodified'] = 1

        # finish
        rows.append(rowdata)

    # turn into data frame and export as CSV
    df = pd.DataFrame(rows) 
    modi_keys = list(modi_keys)
    df.loc[:, modi_keys] = df.loc[:, modi_keys].fillna(0)
    print('\texporting phrase df:', df.shape)
    df.to_csv(paths['phrasedataset'], index=False)

def build_datasets(paths):
    """Load TF and build time / phrase datasets."""

    # load needed TF BHSA data
    TF = Fabric(locations=paths['bhsadata'], silent='deep')
    API = TF.load(
        'kind lex vt pdp ls '
        'lex_utf8 nu code rela '
        'prs prs_gn prs_nu prs_ps '
        'genre mother txt uvf '
    )

    # execute the creation of the data
    time_dataset(paths, API)
    phrase_dataset(paths, API)
