'''
Time objects for analyzing time constructions.
'''

from .helpers import get_index

class Time:
    '''
    The specifications of a Time object
    are stored and made available.
    '''
    
    def __init__(self, cx, tf):
    
        '''
        Tags time constructions
        with specifications found around their 
        time nouns. Returns a dictionary with
        spec strings as keys and 1 as values,
        wherein 1 simply means present.
        '''
        
        # Text-Fabric classes
        F, E, T, L = tf.F, tf.E, tf.T, tf.L
        
        def isQualQuant(word):
            '''
            Determines whether a word is a 
            qualitative quantifier
            '''
            F = tf.F # Text-Fabric Feature Class
            is_quant = F.sem_set.v(word) == 'quant'
            not_card = F.ls.v(word) != 'card'
            if is_quant and not_card:
                return True
            else:
                return False
    
        sent_words = L.d(L.u(cx, 'sentence')[0], 'word')
        times = [time[0] for time in E.role.t(cx) if time[1] == 'time'] or E.nhead.t(L.d(cx, 'phrase')[0])
        quants = [num[0] for num in E.role.t(cx) if num[1] == 'quant']
        adjvs = []
        
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
            m1 = get_index(ph_words, timei-1) # minus 1, etc.
            m2 = get_index(ph_words, timei-2)
            p1 = get_index(ph_words, timei+1) # plus 1, etc.
            p2 = get_index(ph_words, timei+2)
            # relative slots in sentence
            timei_s = sent_words.index(time)
            sp1 = get_index(sent_words, timei_s+1)

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
                            F.lex.v(p1) == 'H' and isQualQuant(p2),
                            isQualQuant(p1) and F.pdp.v(p1) == 'adjv'])
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
            attribs = {'atr', 'adj'}
            subphrases = L.u(time, 'subphrase')
            attr_relas = set(rel_sp for sp in subphrases 
                                 for rel_sp in E.mother.t(sp) 
                                 if F.rela.v(rel_sp) in attribs)

            if attr_relas and  not ({'demon', 'ord', 'attrb', 'qual'} & set(features.keys())):
                features['adjv'] = 1
                adjvs.extend(w for sp in attr_relas for w in L.d(sp, 'word'))

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
        self.adjvs = adjvs
        self.quants = quants or [w for w in L.d(cx, 'word') if isQualQuant(w)]
        self.preps = [ch for ch in L.d(cx, 'chunk') if F.label.v(ch) == 'prep']
        self.cx = cx