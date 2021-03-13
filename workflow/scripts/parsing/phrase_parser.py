import json
from tf.fabric import Fabric
from sly import Lexer, Parser
from sly.lex import Token

class CX:
    """A construction class for representing a word/morpheme."""
    def __init__(self, slot, pos):
        self.slot = slot
        self.cat = pos
    def __repr__(self):
        return f'CX({self.slot})'
        
class BhsaLexer(Lexer):
    """Wrap a BHSA nodes in SLY Token objects.
        
    This enables us to use our custom tokens rather
    than raw string input.
    """

    # initialize standard methods / attributes
    # and add custom ones
    def __init__(self, tf_api, slot2pos):
        super().__init__()
        self.tfa = tf_api
        self.slot2pos = slot2pos
 
    # unique parts of speech values harvested from 
    # the parts of speech parser and pasted here
    # see results/parsing/uniquepos.txt
    tokens = {
         ADJV,
         ADVB,
         ART,
         C,
         CARD,
         CONJ,
         INRG,
#         INTJ,
#         NEGA,
         NOUN,
         ORDN,
         PRDE,
         PREP,
         PRIN,
         PROPN,
         PRPS,
         QUANT,
         SFX,
#         VERB,
    }
    
    def get_token(self, cx, cat, index):
        """Compile custom SLY token object."""
        token = Token()
        token.value = cx
        token.type = cat
        token.index = index
        token.lineno = index
        return token
    
    def tokenize(self, slots):
        """Write over the SLY tokenize method to yield custom data.

        Args:
            slots: iterable of integers which correspond with
                slots / word nodes in BHSA
            slot2pos: a dict that maps slots to their parts of speech
        
        Yields:
            generator of SLY Tokens
        """
        i = 0
        for slot in slots:

            # configure the word
            pos = self.slot2pos[str(slot)]
            cx = CX(slot, pos)
            token = self.get_token(cx, pos, i)
            yield token
            i += 1
            
            # split off morphological tokens
            for morph in self.tokenize_morphs(slot):
                yield self.get_token(morph, morph.cat, i)
                i += 1
        
    def tokenize_morphs(self, slot):
        """Tokenize morphological forms."""
        F = self.tfa.F
        pos = self.slot2pos[str(slot)]
        prs = F.prs.v(slot)
        st = F.st.v(slot)
        if pos != 'PREP':
            if st == 'c':
                yield CX(slot, 'C')
            #elif st == 'a' and prs in {'absent', 'n/a'}:
            #    yield CX(slot, 'A')
        if F.prs.v(slot) not in {'absent', 'n/a'}:
            yield CX(slot, 'SFX')
        
class PhraseParser(Parser):

    """A CFG parser for Biblical Hebrew phrases."""

    # initialize standard methods / attributes
    def __init__(self, error_tracker):
        super().__init__()
        self.error_tracker = error_tracker
    
    tokens = BhsaLexer.tokens

    def error(self, token):
        """Keep track of errors."""
        try:
            self.error_tracker['e'] = token.value.slot
        except:
            self.error_tracker['e'] = None

    precedence = [
        ('right', CONJ),
        ('left', SFX),
    ]
    
    debugfile = 'parser.out'
    
    # -- FINAL PHRASE --
    @_('np', 'pp','cardc',
       'num', 'gp_card', 'gp_num',
       'quant', 'para', 'advb')
    def phrase(self, p):
        if type(p[0]) == int:
            return [p[0]]
        return p[0]

    # -- noun-based phrases --
    @_('NOUN SFX', 'QUANT SFX', 'CARD SFX')
    def np(self, p):
        return p[0].slot
    
    @_('defi', 'appo', 
       'gp', 'adjv')
    def np(self, p):
        return p[0]
    
    # -- definite phrases --
    @_('ART NOUN', 'ART PROPN', 'ART ORDN', 
       'ART PRDE', 'ART QUANT', 'ART CARD',
       'ART ADJV',)
    def defi(self, p):
        return [p[0].slot, p[1].slot, 'DEF']

    @_('ART adjv')
    def defi(self, p):
        return [p[0].slot, p[1], 'DEF']

    # -- prepositional phrases --
    @_('PREP phrase')
    def pp(self, p):
        return [p[0].slot, p[1], 'PP']

    @_('PREP NOUN', 'PREP ADVB', 'PREP PROPN',
       'PREP CARD', 'PREP PRDE', 'PREP QUANT',
       'PREP INRG', 'PREP PRIN', 'PREP PRPS')
    def pp(self, p):
        return [p[0].slot, p[1].slot, 'PP']
    
    @_('PREP SFX')
    def pp(self, p):
        return p[0].slot

    @_('PREP PREP')
    def pp(self, p):
        return [p[0].slot, p[1].slot, 'PP']
    
    # -- genitive phrases --
    @_('NOUN C NOUN', 'NOUN C PROPN',
       'PROPN C NOUN', 'NOUN C QUANT',
       'ADJV C NOUN')
    def gp(self, p):
        return [p[2].slot, p[0].slot, 'GP']
    
    @_('NOUN C np', 'PROPN C np',
       'NOUN C para')
    def gp(self, p):
        return [p[2], p[0].slot, 'GP']

    @_('NOUN C num')
    def gp_num(self, p):
        return [p[2], p[0].slot, 'GP_CARD']
    
    @_('NOUN C CARD')
    def gp_card(self, p):
        return [p[0].slot, p[2].slot, 'GP_NUM']

    # -- quantified phrases -- 
    @_('QUANT C NOUN', 'QUANT C PROPN', 'QUANT C PRDE')
    def quant(self, p):
        return [p[0].slot, p[2].slot, 'QUANT'] 

    @_('QUANT C np')
    def quant(self, p):
        return [p[0].slot, p[2], 'QUANT']
    
    # -- adjectival phrases --
    @_('NOUN ADJV', 'NOUN ADVB', 'NOUN QUANT')
    def adjv(self, p):
        return [p[1].slot, p[0].slot, 'ADJV']

    @_('ADVB NOUN', 'ADVB PRPS')
    def adjv(self, p):
        return [p[0].slot, p[1].slot, 'ADJV']

    @_('ADVB np')
    def adjv(self, p):
        return [p[0].slot, p[1], 'ADJV']

    # -- adverbial phrases -- 
    # !! TO DO: add rule to recognize distributive 
    # constructions suhch as סביב סביב
    @_('ADVB ADVB')
    def advb(self, p):
        return [p[1].slot, p[0].slot, 'ADVB']

    # -- appositional phrases --
    @_('NOUN NOUN', 'NOUN PROPN', 'PROPN PROPN',
       'PROPN NOUN')
    def appo(self, p):
        return [p[1].slot, p[0].slot, 'APPO']
    
    @_('PROPN appo')
    def appo(self, p):
        return [p[1], p[0].slot, 'APPO']

    @_('NOUN gp')
    def appo(self, p):
        return [p[1], p[0].slot, 'APPO']

    @_('defi defi')
    def appo(self, p):
        return [p[1], p[0], 'APPO']

    @_('NOUN defi')
    def appo(self, p):
        return [p[1], p[0].slot, 'APPO']

#    @_('np defi')
#    def appo(self, p):
#        return [p[0], p[0], 'APPO']

    # -- parallel phrases --
    @_('NOUN CONJ NOUN', 'ADVB CONJ ADVB',
       'PROPN CONJ PROPN')
    def para(self, p):
        return [p[0].slot, [p[1].slot, p[2].slot, 'CONJ'], 'PARA']

    @_('NOUN CONJ para', 'PROPN CONJ np', 'NOUN CONJ np')
    def para(self, p):
        return [p[0].slot, [p[1].slot, p[2], 'CONJ'], 'PARA']

    @_('para CONJ NOUN', 'para CONJ PROPN')
    def para(self, p):
        return [p[0], [p[1].slot, p[2].slot, 'CONJ'], 'PARA']

    @_('NOUN C NOUN CONJ gp')
    def para(self, p):
        gp = [p[0].slot, p[2].slot, 'GP']
        cj = [p[3].slot, p[4], 'CONJ']
        return [gp, cj, 'PARA']

    @_('PREP NOUN CONJ pp', 'PREP PROPN CONJ pp')
    def para(self, p):
        pp = [p[0].slot, p[1].slot, 'PP']
        cj = [p[2].slot, p[3], 'CONJ']
        return [pp, cj, 'PARA']

    # -- chained cardinal numbers --
    @_('CARD CARD')
    def cardc(self, p):
        return [p[0].slot, p[1].slot, 'CARDC']
    
    @_('CARD CONJ CARD')
    def cardc(self, p):
        return [p[0].slot, [p[1].slot, p[2].slot, 'CONJ'], 'CARDC']
    
    @_('cardc CONJ cardc')
    def cardc(self, p):
        return [p[0], [p[1].slot, p[2], 'CONJ'], 'CARDC']
    
    @_('CARD CONJ cardc')
    def cardc(self, p):
        return [p[0].slot, [p[1].slot, p[2], 'CONJ'], 'CARDC']
    
    @_('cardc CARD')
    def cardc(self, p):
        return [p[0], p[1].slot, 'CARDC']

    @_('CARD C CARD')
    def cardc(self, p):
        return [p[0].slot, p[2].slot, 'CARDC']
    
    # -- cardinal quantifications --
    @_('NOUN CARD', 'CARD NOUN')
    def num(self, p):
        return [p.CARD.slot, p.NOUN.slot, 'NUM']
    
    @_('cardc NOUN')
    def num(self, p):
        return [p[0], p[1].slot, 'NUM']
    
    @_('cardc np')
    def num(self, p):
        return [p[0], p[1], 'NUM']
    
    @_('CARD np')
    def num(self, p):
        return [p[0].slot, p[1], 'NUM']
    
    @_('np CARD')
    def num(self, p):
        return [p[1].slot, p[0], 'NUM']

    @_('CARD C NOUN')
    def num(self, p):
        return [p[0].slot, p[1].slot, 'NUM']

    @_('CARD C np')
    def num(self, p):
        return [p[0].slot, p[1].slot, 'NUM']

def parse_phrases(samp_path, parsepath, noparsepath, datalocs, API=None):
    """Apply the parser and return metrics on unmatches."""

    # initialize TF
    if API is None:
        TF = Fabric(locations=datalocs['bhsadata'])
        API = TF.load('prs st')

    # load pre-processed parts of speech values
    with open(datalocs['slot2pos'], 'r') as infile:
        slot2pos = json.load(infile)

    # load samples
    with open(samp_path, 'r') as infile:
        samples = json.load(infile)

    # initialize lexer/parser
    lexer = BhsaLexer(API, slot2pos)

    # initialize parser
    # error_tracker allows us to save
    # an offending value and keep it associated with 
    # the respective phrase node number;
    # gets reset during each iteration of the loop
    error_tracker = {'e': None}
    parser = PhraseParser(error_tracker)

    # run the parser
    parsed = {}
    errors = {}
    for ph_node, slotset in samples.items():

        # do not run parser for single-word phrases
        if len(slotset) == 1:
            parsed[ph_node] = slotset
            continue

        # parse multi-word phrases
        tokens = list(lexer.tokenize(slotset))
        parsing = parser.parse(t for t in tokens)
        if parsing is not None:
            parsed[ph_node] = parsing
        else:
            toks = [(t.type, t.value.slot) for t in tokens]
            errors[ph_node] = (error_tracker['e'], str(toks))

            # reset error tracker
            error_tracker['e'] = None

    # export
    with open(parsepath, 'w') as outfile:
        json.dump(parsed, outfile, indent=2)

    with open(noparsepath, 'w') as outfile:
        json.dump(errors, outfile, indent=2)
