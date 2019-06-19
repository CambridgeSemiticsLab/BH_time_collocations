'''
Time objects for analyzing time constructions.
'''

from __main__ import F, E, T, L

def getIndex(thislist, index):
    '''
    A safe way to get index from 
    a list/tuple. If indexError returns None.
    '''
    try:
        return thislist[index]
    except IndexError:
        return None
    
def isQualQuant(word):
    if F.sem_set.v(word) == 'quant' and F.ls.v(word) != 'card':
        return True
    else:
        return False
    
class Time:
    '''
    The specifications of a Time object
    are stored and made available.
    '''
    
    def __init__(self, cx):
    
        '''
        Tags time constructions
        with specifications found around their 
        time nouns. Returns a dictionary with
        spec strings as keys and 1 as values,
        wherein 1 simply means present.
        '''
    
        sent_words = L.d(L.u(cx, 'sentence')[0], 'word')
        times = [time[0] for time in E.role.t(cx) if time[1] == 'time'] or E.nhead.t(L.d(cx, 'phrase')[0])
        quants = [num[0] for num in E.role.t(cx) if num[1] == 'quant']
        
        features = {}

        for time in times:
            
            phrase = L.u(time, 'phrase')[0]
            ph_words = L.d(phrase, 'word')
            dep_cl = next((cl for cl in E.mother.t(phrase) if F.rela.v(cl) == 'Attr'), None)
            
            # phrase type, PP or NP, wherein AdvP are considered a kind of NP
            typ = 'PPtime' if F.typ.v(phrase) == 'PP' else 'time'
            features[typ] = 1

            # get relative slot positions
            timei = ph_words.index(time)
            m1 = getIndex(ph_words, timei-1) # minus 1, etc.
            m2 = getIndex(ph_words, timei-2)
            p1 = getIndex(ph_words, timei+1) # plus 1, etc.
            p2 = getIndex(ph_words, timei+2)
            # relative slots in sentence
            timei_s = sent_words.index(time)
            sp1 = getIndex(sent_words, timei_s+1)

            # preceding article
            if F.lex.v(m1) == 'H':
                features['H'] = 1

            # plurals
            if F.nu.v(time) == 'pl' and F.pdp.v(time) != 'prde':
                features['pl'] = 1

            elif F.nu.v(time) == 'du':
                features['quant'] = 1
                features['du'] = 1

            # pronom suffixs
            if F.prs.v(time) not in {'absent', 'n/a'}:
                features['sffx'] = 1

            # check quant & qual quants
            if quants:
                features['quant'] = 1
                features['card'] = 1

            is_qualq = any([isQualQuant(m1),
                            isQualQuant(m2) and F.lex.v(m1) == 'H', # attr patterns
                            F.lex.v(p1) == 'H' and isQualQuant(p2)])
            if is_qualq:
                features['quant'] = 1
                features['qual'] = 1

            # constructs
            if F.st.v(time) == 'c':
                next_word = p1 if F.pdp.v(p1) != 'art' else p2
                features['construct'] = 1
                if F.pdp.v(sp1) == 'verb':   
                    features['cons+VC'] = 1
                elif F.pdp.v(next_word):
                    features['cons+NP'] = 1

            # h.time.h.spec pattern
            if F.lex.v(m1) == 'H' and F.lex.v(p1) == 'H':
                features['attr_patt'] = 1

            # demonstrative / ordinal / qualquant / spec
            is_dem = any([F.pdp.v(p1) == 'prde',
                          F.lex.v(p1) == 'H' and F.pdp.v(p2) == 'prde'])
            is_ordn = F.lex.v(p1) == 'H' and F.ls.v(p2) == 'ordn'

            is_spec = all([F.lex.v(m1) == 'H', 
                           F.lex.v(p1) == 'H',
                           not F.lex.v(p1) == 'H' and isQualQuant(p2),
                           not is_dem, not is_ordn, not is_qualq])
            if is_dem:
                features['demon'] = 1
            elif is_ordn:
                features['ord'] = 1
            elif is_spec:
                features['attrb'] = 1

            # attributives
            small_sp = next(iter(sorted(L.u(time, 'subphrase'))), 0)
            attr_relas = set(rel_sp for rel_sp in E.mother.t(small_sp) if F.rela.v(rel_sp) in {'atr', 'adj'})

            if attr_relas and  not {'demonstrative', 'ordinal', 'attribute', 'qualitative quant.'} & set(features.keys()):
                features['adjv'] = 1

        # tag relative/attributive specs dependent on phrase
        if dep_cl:
            rel = 'rela' if 'Rela' in set(F.function.v(ph) for ph in L.d(dep_cl, 'phrase')) else ''
            clkind = F.kind.v(dep_cl)        
            features[f'{rel}+{clkind}'] = 1        

        tag = '.'.join(features.keys())
        result = (cx,) + L.d(cx, 'phrase') + tuple(times)
    
        self.tag = tag
        self.result = result 
        self.specs = features
        self.times = times
        self.quants = quants or [w for w in L.d(cx, 'word') if isQualQuant(w)]
        self.cx = cx