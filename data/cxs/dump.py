# load construction objects from TF resource

import sys
import pickle
import collections
from datetime import datetime
import paths 
from word_grammar import Words
from phrase_grammar import Subphrases, Phrases
from phrase_classes import SinglePhrase
from dataset import build_dataset
from tf_tools.load import load_tf

# load semantic vectors
with open(paths.semvector, 'rb') as infile: 
    semdist = pickle.load(infile)

def load_cxs(tf_api, semdist, debug=False):
    """Load constructional objects"""
    
    A = tf_api
    A.api.makeAvailableIn(globals())    

    # retrieve phrases
    def disjoint(ph):
        """Isolate phrases with gaps."""
        ph = L.d(ph,'word')
        for w in ph:
            if ph[-1] == w:
                break
            elif (ph[ph.index(w)+1] - w) > 1:
                return True

    alltimes = [
        ph for ph in F.otype.s('timephrase') 
    ]
    
    timephrases = [ph for ph in alltimes if not disjoint(ph)]

    print(f'{len(timephrases)} phrases ready')


    # -- WORDS -- 
    
    words = Words(A) # word CX builder

    # analyze all matches; return as dict
    start = datetime.now()
    print()
    print(f'{datetime.now()-start} beginning word construction analysis...')
    wordcxs = words.cxdict(
        s for tp in timephrases
            for s in L.d(tp,'word')
    )
    print(f'\t{datetime.now() - start} COMPLETE \t[ {len(wordcxs)} ] words loaded')
    
    # -- SUBPHRASES -- 
    
    # time phrase CX builder
    spc = Subphrases(wordcxs, A, semdist)
    
    phrase2cxs = collections.defaultdict(list)
    nocxs = []

    # time it
    print()
    start = datetime.now()
    print(f'{datetime.now()-start} beginning subphrase analysis...')

    for i, phrase in enumerate(timephrases):

        # analyze all known relas
        elements = L.d(phrase,'word')

        # analyze with debug exceptions
        try:
            cxs = spc.analyzestretch(elements)
        except:
            sys.stderr.write(f'\nFAIL...running with debug...\n')
            pretty(phrase)
            spc.analyzestretch(elements, debug=debug)
            raise Exception('...debug complete...')

        # save those phrases that have no matching constructions
        if not cxs:
            nocxs.append(phrase)
        else:
            phrase2cxs[phrase] = cxs

        # report status
        if i % 500 == 0 and i:
            print(f'\t{datetime.now()-start}\tdone with iter {i}/{len(timephrases)}')

    print(f'{datetime.now()-start}\tCOMPLETE')
    print('-'*20)
    print(f'{len(phrase2cxs)} phrases matched with Constructions...')
    print(f'{len(nocxs)} phrases not yet matched with Constructions...')
    

    print('Classifying single timephrases...')

    print('building preprocess data...')
    # compile acceptable head lexemes from single-phrased CXs
    good_heads = set()
    for ph, cx_data in phrase2cxs.items():
        if len(cx_data) == 1:
            cx = cx_data[0]
            head = list(cx.getsuccroles('head'))[-1]
            good_heads.add(F.lex.v(head))
    
    # tag the time cxs with classifications
    sp = SinglePhrase(phrase2cxs.values(), good_heads, A)  
    sp.label_cxs()

    return {
        'wordcxs': wordcxs, 
        'phrase2cxs': phrase2cxs,
        'class2cx': sp.class2cx,
    } 

# load up Text-Fabric
TF, api, A = load_tf()
print()

if '-lite' not in sys.argv:
    # -- Run Construction Analysis -- 
    cxs = load_cxs(A, semdist)

    # -- Dump Construction Objects --
    with open(paths.cxs, 'wb') as outfile:
        pickle.dump(cxs, outfile)
    print()
    print(f'Dumping cxs into {paths.cxs}')

# read back in for primary data build
from cx_analysis.load import cxs

# -- Compile Primary CSV Dataset -- 
print(f'Building primary dataset')
constructions = [c for cs in cxs['phrase2cxs'].values() for c in cs]
dataset = build_dataset(constructions, A.api)
dataset.to_csv(paths.main_table, sep='\t')
print(f'DONE! Pushed to {paths.main_table}')
