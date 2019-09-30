class Nominals:
    '''
    Prepares a set of words which are 
    neither verbs nor articles nor conjunctions.
    I.e. this incldudes "nouns" in the traditional 
    sense, but also "adjectives" and "adverbs".
    It also includes participles.
    '''
    
    def __init__(self, tf, **wsets):
        
        F = tf.api.F
        preps = wsets['preps']
        
        sps = {'subs', 'nmpr', 'adjv', 'advb'}
        
        self.nominals = set()
        
        for w in F.otype.s('word'):
            
            # skip prepositional words
            if w in preps:
                continue
                
            # add nominals
            if F.sp.v(w) in sps:
                self.nominals.add(w)
            elif (F.sp.v(w) == 'verb' and F.vt.v(w) in {'ptcp', 'ptca'}):
                self.nominals.add(w)