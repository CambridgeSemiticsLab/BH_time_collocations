import numpy as np
import pandas as pd
from tf_tools.load import load_tf
import tf_tools.formatting as form

print("Beginning export of broad sample...\n")

# set up Text-Fabric methods
TF, API, A = load_tf('nhead')
F, E, T, L = A.api.F, A.api.E, A.api.T, A.api.L

# gather the dataset
phrase_dataset = []

for phrase in F.otype.s('phrase'):
    
    # exclude Aramaic portions
    lang = F.language.v(L.d(phrase,'word')[0])
    if lang != 'Hebrew':
        continue
    
    book, chapter, verse = T.sectionFromNode(phrase)
    book_sbl = form.book2sbl[book]
    sentence = L.u(phrase, 'sentence')[0]
    ref = f'{form.book2sbl[book]} {chapter}:{verse}'
    function = F.function2.v(phrase)
    s_function = form.simplified_functions.get(function, function)
    n_words = len(L.d(phrase, 'word'))
    n_phrase_atoms = len(L.d(phrase, 'phrase_atom'))
    heads = E.nhead.t(phrase)
    n_heads = len(heads)
    first_head = heads[0]
    head_lexs = '|'.join(F.lex.v(h) for h in heads)
    head_utf8 = '|'.join(T.text(h).strip() for h in heads)
    head_lexnodes = tuple(L.u(h,'lex')[0] for h in heads)
    daughters = E.mother.t(phrase)
    n_daughters = len(daughters)
    mothers = E.mother.f(phrase)
    n_mothers = len(mothers)
    d_relas = '|'.join(F.rela.v(d) for d in daughters)
    ph_typ = F.typ.v(phrase)
    rela = F.rela.v(phrase)
    time_phrase = L.u(phrase, 'timephrase')
    
    phrase_dataset.append({
        'node': phrase,
        'ref': ref,
        'book': book_sbl,
        'text': T.text(phrase),
        'sentence': T.text(sentence),
        'type': ph_typ,
        'function2': function,
        's_function': s_function,
        'rela': rela,
        'n_words': n_words,
        'n_phrase_atoms': n_phrase_atoms,
        'n_heads': n_heads,
        'head_lex': head_lexs,
        'head_lex_nodes': head_lexnodes,
        'n_daughters': n_daughters,
        'daught_relas': d_relas or np.nan,
        'n_mothers': n_mothers,
        'n_relas': n_daughters + n_mothers,
        'in_timephrase': bool(time_phrase),
    })

    if len(phrase_dataset) % 10000 == 0:
        print(f"\t{len(phrase_dataset)} phrases analyzed...")

print("exporting...")

phrase_df = pd.DataFrame(phrase_dataset)
phrase_df.set_index('node', inplace=True)
phrase_df.to_csv('broad_dataset.csv') 

print('DONE!')
