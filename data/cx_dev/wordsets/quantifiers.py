from positions import Positions

class Quants:
    
    def __init__(self, tf):
        '''
        Input: instance of Text Fabric
        Output: a series of sets containing
        custom quantifier nodes
        '''

        F = tf.api.F
        
        quants = set()
        quants |= set(F.ls.s('card')) & set(F.otype.s('word'))

        # -- contextual cases -- 
        
        for w in F.otype.s('word'):
        
            P = Positions(w, 'phrase_atom', tf).get
        
            # -- custom lexemes-- 
            custom = {'KL/', 'M<V/', 'JTR/',
                      'M<FR/', 'XYJ/', '<FRWN/',
                      'C>R=/', 'MSPR/', 'RB/', 
                      'RB=/', 'XMJCJT/'}
            if F.lex.v(w) in custom:
                quants.add(w)
                continue

            # -- the Hebrew idiom: בנ + quantifier for age --
            conds = [
                P(0, 'lex') == 'BN/',
                P(0, 'st') == 'c',
                P(0, 'nu') == 'sg',
                P(1, 'ls') == 'card',
            ]
            if all(conds):
                quants.add(w)
                continue
        
        self.quants = quants