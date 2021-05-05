import pandas as pd
from tf.fabric import Fabric
from tools.load_parse import ParseLoader
from .tokenizers import tokenize_lexemes
from .nav_time import get_times
from .verb_form import get_verbform

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
        'lex_utf8'
    )
    F, E, T, L = API.F, API.E, API.T, API.L

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
        cl_text = T.text(clause)
        sentn = L.u(clause, 'sentence')[0]
        sent_text = T.text(sentn)
        bk, ch, vs = T.sectionFromNode(versen)
        verse = f'{bk} {ch}:{vs}'
        name = data['functions'][0]
        function = function_simp.get(name, name)
        quality = quality_map.get(function, None)
        times = sorted(get_times(data))
        time_lexs = '|'.join(F.lex.v(t) for t in times)
        lex_token = tokenize_lexemes(
            ph_parses, 
            API,
            heads=times,
        )
        if len(times) == 1:
            is_advb = F.pdp.v(times[0]) == 'advb'
        else:
            is_advb = False        
  
        data = {
            'node': clause,
            'book': bk,
            'verse': verse,
            'function': function,
            'quality': quality,
            'name': name,
            'text': text,
            'time_lexs': time_lexs,
            'n_times': len(times),
            'lex_token': lex_token,
            'is_advb': is_advb,
            'cl_kind': F.kind.v(clause),
            'clause': cl_text,
            'sentence': sent_text,
        }
        data.update(
            get_clause_data(clause, API)
        )
        rows.append(data)

    # turn into data frame and export as CSV
    df = pd.DataFrame(rows) 
    print('exporting:', df.shape)
    df.to_csv(paths['dataset'], index=False)
