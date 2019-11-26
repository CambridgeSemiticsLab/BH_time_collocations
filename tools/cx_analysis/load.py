# load construction objects from TF resource

import collections
from datetime import datetime
from locations import semvector
from .grammar import Words
from .phrase_grammar import Subphrases, Phrases

# load semantic vectors
with open(semvector, 'rb') as infile: 
    semdist = pickle.load(infile)

def load_cxs(tf_api, debug=False):
    """Load constructional objects"""
    
    A = tf_api
    
    # -- WORDS -- 
    
    words = Words(A) # word CX builder

    # analyze all matches; return as dict
    start = datetime.now()
    print(f'Beginning word construction analysis...')
    wordcxs = words.cxdict(
        s for tp in timephrases
            for s in L.d(tp,'word')
    )
    print(f'\t{datetime.now() - start} COMPLETE \t[ {len(wordcxs)} ] words loaded')
    
    # -- SUBPHRASES -- 
    
    # time phrase CX builder
    spc = Subphrases(wordcxs, semdist, A)
    
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
    
    