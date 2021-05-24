import json
import pickle
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
from .tag_args import clause_args, tag_position
from .clause_tree import get_successors

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

def get_clause_data(clause, API, 
                    parsedata, 
                    do_args=True,
                    count_succs=True,
                    ):
    """Build data on a clause."""

    # TF methods
    F, E, T, L = API.F, API.E, API.T, API.L

    cl_data = {}

    # collect data of clause, but also
    # of embedding linguistic unit like sentence
    firstw = L.d(clause, 'word')[0]
    firstw_fmts = get_word_formats([firstw], parsedata['slot2pos'], API)
    cl_data['firstw'] = firstw_fmts['latex']
    versen = L.u(firstw, 'verse')[0]
    genre = cl_data['genre'] = F.genre.v(versen)
    domain = cl_data['domain'] = permissive_q(clause, API)

    # process genre/domain to main genres
    main_gens = {'prose', 'prophetic', 'poetry', 'instruction'}
    main_doms = {'Q', 'N'}
    if genre in main_gens and domain in main_doms:
        if genre == 'prophetic':
            genre = 'prophecy'
        if genre == 'prose':
            cl_data['main_genre'] = f'{genre}-{domain}'
        else:
            cl_data['main_genre'] = genre

    cl_data['gendom'] = cl_data['genre'] + '.' + cl_data['domain']
    cl_data['clause'] = T.text(clause)
    sentn = L.u(clause, 'sentence')[0]
    cl_data['sentence'] = T.text(sentn)
    cl_data['cl_rela'] = clause_relator(clause, API)
    cl_data['cl_type'] = F.typ.v(clause)
    cl_data['cl_kind'] = F.kind.v(clause)

    # build data for clause verbs
    verbs = [
        w for w in L.d(clause,'word')
            if F.pdp.v(w) == 'verb'
    ]
    verb = verbs[0] if verbs else None
    cl_data['verb'] = verb
    cl_data['verb_etcbc'] = F.lex.v(verb) if verb else None
    cl_data['verb_utf8'] = F.lex_utf8.v(verb)
    verb_text = T.text(verb, fmt='text-orig-plain') if verb else None
    if verb:
        cl_data['verb_text'] = verb_text.strip()
        cl_data['verbform'] = get_verbform(
            verb, 
            API,
            parsedata['bhsa2gbi'] 
        )
        cl_data['verb_stem'] = F.vs.v(verb)
    else:
        cl_data['verbform'] = 'Ø'

    # build data for clause arguments
    if do_args:
        refword = verb if verb else firstw 
        cl_args = clause_args(
            refword, 
            API,
            parsedata['functions'],
        )
        cl_data['cl_args'] = cl_args
        cl_data['has_objc'] = 1*('O' in cl_args)
        cl_data['has_cmpl'] = 1*('C' in cl_args)
        cl_data['has_subj'] = 1*('C' in cl_args)
        cl_data['has_oc'] = 1*(cl_data['has_objc'] or cl_data['has_cmpl'])
        cl_data['Time Position'] = tag_position(cl_args)

    # add-on cl type that accounts for wayehi / wehaya
    cl_type2 = cl_data['cl_type']
    verb_form = cl_data.get('verbform')
    verb_text = cl_data.get('verb_text')
    if verb_form == 'wayq' and verb_text == 'יהי':
        cl_type2 = 'WayH'
    elif verb_form == 'wqtl' and verb_text == 'היה':
        cl_type2 = 'WQtH'
    cl_data['cl_type2'] = cl_type2

    # count successors based on the first clause atom
    if count_succs:
        clatom = L.d(clause, 'clause_atom')[0]
        succs = list(get_successors(clatom, API))
        # filter to main clause succs
        main_cls = set(
            L.u(ca, 'clause')[0] for ca in succs
        )
        main_cls = [c for c in main_cls if c != clause]
        main_succs = [
            c for c in main_cls
                if clause_relator(c, API) == 'Main'
        ]
        cl_data['cl_nsuccs'] = len(main_succs)

    # return clause data
    return cl_data

def get_word_formats(words,
                     slot2pos, 
                     API, 
                     prefix='', joiner='|', 
                     latexrow='latex'):
    """Retrieve word formats for a list of words."""

    # BHSA methods
    F = API.F

    # format data on heads
    words_etcbc = f'{joiner}'.join(F.lex.v(w) for w in words)
    words_utf8 = f'{joiner}'.join(F.lex_utf8.v(w) for w in words)
    words_utf8d = (
        get_display(remove_shindots(words_utf8))
    )
    # formatting for inclusion in latex docs (esp. tables)
    words_latex = f'{joiner}'.join(
        '\texthebrew{%s}' % F.lex_utf8.v(w)
            for w in reversed(words)
    )
    words_latex = remove_shindots(words_latex)

    # formatting using custom-made values
    words_pos = f'{joiner}'.join(
        slot2pos[w] for w in words
    )

    # package and ship
    formats = {
        f'{prefix}etcbc': words_etcbc,
        f'{prefix}utf8': words_utf8,
        f'{prefix}utf8d': words_utf8d,
        f'{prefix}POS': words_pos,
        f'{latexrow}': words_latex,
    }
    return formats

def add_tense(timedata, modifiers):
    """Add tense designations based on several rules."""

    # ignore if tense already recorded
    if timedata.get('tense'):
        return {}
    
    # calculate present tense constructions
    pres_rules = (
        timedata['front'] == 'Ø'
        and bool(modifiers.get('DEF'))
        and len(modifiers) == 2
    )
    if pres_rules:
        return {'tense': 'PRES'}

    # calculate future tense constructions
    fut_rules = (
        timedata['function'] == 'simultaneous'
        and timedata.get('demon_type') == 'THAT'
        and timedata.get('verbtense') in {'FUT', 'MOD shall'}
    )
    if fut_rules:
        return {'tense': 'FUT'}

    # calculate past tense constructions
    past_rules = (
        timedata['function'] == 'simultaneous'
        and timedata.get('demon_type') == 'THAT'
        and timedata.get('verbtense') in {'PAST'}
    )
    if past_rules:
        return {'tense': 'PAST'}

    # at this point there is no matches; return empty dict
    return {} 

def time_dataset(paths, parsedata, API):
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
    time_data = parsedata['times']
    phrase_data = parsedata['phrases']

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
                parsedata['slot2pos'],
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
        tense = data.get('tenses', [[None, None]])[0][0]
        lex_token = tokenize_lexemes(
            ph_parses, 
            API,
            heads=times,
        )

        # get linguistic head
        first_parse = ph_parses[0]
        rowdata['head_utf8'] = F.lex_utf8.v(nt.get_head(first_parse))

        notadvbs = {
            'LJLH/', 'BQR=/', 'JWM/'
        }
        if (len(times) == 1) and rowdata['times_etcbc'] not in notadvbs:
            is_advb = 1*(F.pdp.v(times[0]) == 'advb')
        else:
            is_advb = 0 

        main_functions = {
            'simultaneous',
            'atelic_ext',
            'telic_ext',
            'anterior_dur',
            'posterior',
            'posterior_dur',
            'habitual',
            'anterior',
            'dist_fut',
            'dist_past',
        }
        rowdata['mainfunction'] = 1*(function in main_functions)

        rowdata.update({
            'function': function,
            'quality': quality,
            'name': name,
            'text': text,
            'n_times': len(times),
            'lex_token': lex_token,
            'is_advb': is_advb,
            'tense': tense,
        })

        # add clause-based data
        rowdata.update(
            get_clause_data(clause, API, parsedata)
        )

        rowdata['cl_clust50'] = (
            parsedata['clclusters']['50']['clusters'][str(clause)]
        )
        rowdata['cl_clust10'] = (
            parsedata['clclusters']['10']['clusters'][str(clause)]
        )

        # get various additional data centered on the verb
        if rowdata.get('verb'):
            verb = rowdata['verb']
            verbtense = parsedata['tenses'].get(verb,{}).get('esv_TAMsimp')
            
            # imperatives can be missed by the English alignment+parsing process;
            # but semantically these are straightforward; we fill in the
            # gap when this happens by simply assigning the sense as imperative
            # as it seems unlikely the English would often choose another sense
            if rowdata.get('verbform') == 'impv' and not verbtense:
                verbtense = 'IMPV'

            # add the tense
            rowdata['verbtense'] = verbtense
        
            # calc clause args
            cl_args = rowdata['cl_args']
            rowdata['vt_order']  = ''.join(
                l for l in cl_args if l 
                    if l in {'T', 'V'}
            )
       
        # add modifier data
        # this data will typically be used to analyze single-phrased TAs
        # so we just take the first ph_parse
        modifiers = get_modis(
            ph_parses[0], 
            API, 
            boolean=False
        )
        modi_keys |= set(modifiers)
        bool_mods = {m:1 for m in modifiers}
        rowdata.update(
            bool_mods
        )

        # process phrase types
        if modifiers.get('ØPP') and is_advb:
            rowdata['front'] = rowdata['times_utf8']
            rowdata['ph_type'] = 'ADVB'
        elif modifiers.get('ØPP'):
            rowdata['front'] = 'Ø'
            rowdata['ph_type'] = 'NP'
        elif modifiers.get('PP'):
            pp_strings = get_word_formats(
                sorted(modifiers['PP']),
                parsedata['slot2pos'],
                API,
                joiner='+'
            )
            rowdata['front'] = pp_strings['utf8']
            rowdata['ph_type'] = 'PP'

        # build semantic/formal tag
        qual_map = {
            'location': 'l',
            'duration': 'd',
            'sequence': 's',
            'iteration': 'i',
        }
        no_lex = {'LJLH/', 'JWMM', 'JWM/'}
        if rowdata.get('front_etcbc') and rowdata.get('quality'):
            front = rowdata['front_etcbc']
            qual = rowdata['quality'] 
            if front == 'advb':
                front = rowdata['times_etcbc']
            front = front.replace('/', '').replace('=', '')
            qual_tag = qual_map[qual]
            rowdata['tag'] = f'{qual_tag}_{front}'

        # get demonstrative data
        demon_map = {
            "Z>T": "THIS",
            "HJ>": "THAT",
            "HMH": "THAT",
            ">LH": "THIS",
            "HM": "THAT",
            "HW>": "THAT",
            "ZH": "THIS",
        }    
        if demon := modifiers.get('DEMON'):
            rowdata['demon_type'] = demon_map[F.lex.v(demon[0])]
   
        # mark unmodified words as adverbs
        if not modifiers or (len(modifiers) == 1 and {'PP', 'ØPP'} & set(modifiers)):
            rowdata['unmodified'] = 1
        else:
            rowdata['unmodified'] = 0

        # add tense data based on modifiers and other feats
        rowdata.update(
            add_tense(rowdata, modifiers)
        )

        rowdata['has_time'] = 1

        # finish
        rows.append(rowdata)

    # turn into data frame and export as CSV
    df = pd.DataFrame(rows) 
    modi_keys = list(modi_keys)
    df.loc[:, modi_keys] = df.loc[:, modi_keys].fillna(0)
    print('\texporting time df:', df.shape)
    df.to_csv(paths['timedataset'], index=False)

def phrase_dataset(paths, parsedata, API):
    """Construct tabular dataset of time adverbials."""

    print('Building phrase dataset...')

    # text-fabric methods
    F, E, T, L = API.F, API.E, API.T, API.L

    # load data parsed in this project
    phrase_data = parsedata['phrases']
    functions = parsedata['functions']

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
                parsedata['slot2pos'],
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
        phatoms = [
            pa for ph in phrases
                for pa in L.d(ph,'phrase_atom')
        ]
        head_ph = L.u(heads[0], 'phrase')[0]
        text = T.text(slots)
        function = functions[head_ph]
        word_lexs = ' '.join(
            F.lex.v(s) for s in slots
        )
        typs = '|'.join(F.typ.v(pa) for pa in phatoms)
        rowdata.update({
            'function': function,
            'text': text,
            'types': typs,
            'n_heads': len(heads),
            'word_lexs': word_lexs,
            'n_words': len(slots),
            'n_phatoms': len(phatoms),
        })

        # add clause-based data
        clause = L.u(head_ph, 'clause')[0]
        rowdata['clause_node'] = clause
        rowdata.update(
            get_clause_data(
                clause, API, 
                parsedata, 
                do_args=False,
                count_succs=False,
            )
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

def clause_dataset(paths, parsedata, API):
    """Build data on all clauses that are not already parsed."""
    F, L = API.F, API.L
    
    print('Building clause dataset...')

    rows = []

    for clause in F.otype.s('clause'):
    
        # skip clauses that are already parsed
        if clause in parsedata['times']:
            continue

        row_data = {'node': clause, 'tag': 'NT', 'has_time': 0} # no time clauses
        row_data.update(get_clause_data(clause, API, parsedata))
        rows.append(row_data)

    df = pd.DataFrame(rows)
    print('\texporting clause df:', df.shape)
    df.to_csv(paths['clausedataset'], index=False)

def build_datasets(paths):
    """Load TF and build time / phrase datasets."""

    # load various data for feature processing
    data = {
        'slot2pos': ParseLoader(paths['slot2pos']).load(),
        'times': ParseLoader(paths['timedata']).load(),
        'tenses': ParseLoader(paths['tensedata']).load(),
        'phrases': ParseLoader(paths['phrasedata']).load(),
        'functions': ParseLoader(paths['functions']).load(),
    }

    with open(paths['clclusters'], 'r') as infile:
        data['clclusters'] = json.load(infile)

    with open(paths['bhsa2gbi_verb'], 'rb') as infile:
        data['bhsa2gbi'] = pickle.load(infile)
    
    # load needed TF BHSA data
    TF = Fabric(locations=paths['bhsadata'], silent='deep')
    API = TF.load(
        'kind lex vt pdp ls '
        'lex_utf8 nu code rela '
        'prs prs_gn prs_nu prs_ps '
        'genre mother txt uvf typ '
        'g_prs_utf8 sp vs typ '
    )

    # execute the creation of the data
    time_dataset(paths, data, API)
    phrase_dataset(paths, data, API)
    clause_dataset(paths, data, API)
