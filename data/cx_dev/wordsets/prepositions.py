'''
Delivers custom prepositions.
'''

import csv
import collections
from positions import Positions

class Preps:
    
    def __init__(self, tf):
        
        F, L = tf.api.F, tf.api.L
        self.preps = set(w for w in F.otype.s('word') if F.pdp.v(w)=='prep')        
        
        # -- ETCBC "potential preps" --
        # דרך is excluded since this is speculative
        self.preps |= set(w for w in F.otype.s('word')
                             if F.ls.v(w) == 'ppre'
                             and F.lex.v(w) != 'DRK/'
                         )
    
        
        # -- context-specific prepositions --
        
        # prepare manually annotated contexts
        # - ראשׁ
        with open('../annotations/rosh_constructs.csv', 'r') as infile:
            rosh_annotes = list(csv.reader(infile))[1:]
            rosh_kids = [row[-2] for row in rosh_annotes if row[-1] == 'y']
        
        # begin word loop
        for w in F.otype.s('word'):
            
            # skip unnecessary pos
            if F.sp.v(w) not in {'prep', 'subs', 'adjv', 'advb'}:
                continue
            
            P = Positions(w, 'phrase_atom', tf).get
            
            # Set up potential construct position
            # NB! still requires a check for construct state
            construct = (
                P(1)
                    if P(1, 'sp') != 'art'
                    else None
                or P(2)
                    if P(1, 'sp') == 'art'
                    else None
            )
            
            # -- lexemes that must be preceded by another prep -- 
            # These sets could benefit from further investigation
            preprep = {'PNH/', 'TWK/', 'QY/', 'QYH=/', 'QYT/', '<WD/'}
            conds = [
                
                P(0, 'lex') in preprep,
                P(0, 'prs') == 'absent',
                (P(-1, 'pdp') == 'prep'
                     or P(-1, 'ls') == 'ppre'
                     or P(-1, 'lex') == 'KL/'),
                P(1, 'lex') != 'W'
                
            ]
            if all(conds):
                self.preps.add(w)
                continue
                
            # -- בד "alone" when preceded by ל, with a meaning of "except" --
            elif (F.lex.v(w) == 'BD/' and P(-1, 'lex') == 'L'):
                self.preps.add(w)
                continue

            # -- cases of אחרית -- 
            '''
            Some are prep and some are substantival.
            e.g. subs. אחרית רשׁעים "end of evil doers" (Ps 37:38)
            others are used prepositionally to indicate position
            the semantics of the components are important for determining 
            which sense is employed. 
            All cases in Time function phrases appear prepositional.
            Cases with animate nouns appear to be substantival;
            those cases can be manually excluded with a lexeme exclusion.
            NB: גים in Jer 50:12 is used non-personally and thus not excluded
            Excluded: איוב and רשעים
            '''
            conds = [
                P(0, 'lex') == '>XRJT/',
                P(0, 'lex') == 'c',
                F.lex.v(construct) not in {'>JWB/', 'RC</'}
            ]
            if all(conds):
                self.preps.add(w)
                continue
                
            # -- prepositional ראשׁ --        
            '''
            רֹאשׁ "head" is a word seems to become prepositional when attached to 
            inanimate nouns. However, BHSA currently has no way of separating animate 
            and inanimate nouns. This is a case where my current research on noun embeddings could 
            soon provide the necessary data. The method applied there distinguishes animate from inanimate 
            nouns with fairly good success. Eventually it could be possible to identify prepositions 
            like ראשׁ by searching for words in construct to it with an animacy rating below a 
            given threshold. However, since the embeddings do not yet account for various word senses, 
            more work is needed before more reliable results can be procured. 
            For now, this case is tagged as one which in the near future might 
            be improved with empirical semantic data.

            The possibility of automatically distinguishing prepositional cases of ראשׁ also opens the 
            possibility for other similar body part prepositional uses such as פֶּה "mouth", יָד "hand", שָׂפָה "lip", 
            and potentially others.

            For now, a brute force method is applied which filters cases out by testing validity of the lexeme 
            subservient to ראשׁ. The lexemes have been manually selected by examining a set of noun lexemes in construct 
            to ראשׁ and by double checking the context. The process is accomplished through a simple manual annotation. 
            A spreadsheet has been exported and hand annotated. The results are drawn in and applied to ראשׁ.
            '''
            conds = [
                P(0, 'lex') == 'R>C/',
                P(0, 'st') == 'c',
                P(-1, 'pdp') == 'prep',
                (F.lex.v(construct) in rosh_kids
                    or P(1, 'ls') == 'card')
            ]
            if all(conds):
                self.preps.add(w)
                continue
                
            # ETSEM (time phrase)
            
            conds = [
                P(0,'lex') == '<YM/',
                P(0,'st') == 'c',
                F.function.v(L.u(w,'phrase')[0]) == 'Time'
            ]
            
            if all(conds):
                self.preps.add(w)
                continue