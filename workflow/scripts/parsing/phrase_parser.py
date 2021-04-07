import json
from tf.fabric import Fabric
from sly import Lexer, Parser
from sly.lex import Token
from tools.positions import PositionsTF

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
    def __init__(self, tf_api, slot2pos, word_seps):
        super().__init__()
        self.tf_api = tf_api
        self.slot2pos = slot2pos
        self.word_seps = word_seps
 
    # unique parts of speech values harvested from 
    # the parts of speech parser and pasted here
    # see results/parsing/uniquepos.txt
    tokens = {
         ADJV,
         ADVB,
         ART,
         CARD,
         CARD1,
         CONJ,
         CONJPP,
         CONJGP,
         CONJADJV,
         CONJCARD,
         CONJQUANT,
         INRG,
         NOUN,
         ORDN,
         PRDE,
         PREP,
         ADJVPREP,
         PRIN,
         PROPN,
         PRPS,
         QUANT,
         NEGA,
         GAM,
         C, # construct state
         A, # absolute state
         SFX, # pronominal suffixes
         SEP, # separator for [np np] disambiguation
         SEPPREP,
         SEPGP,
    }
    
    def getP(self, node, context='phrase_atom'):
        """Get Positions object for a TF node."""
        return PositionsTF(node, context, self.tf_api).get

    def get_token(self, cx, cat, index):
        """Compile custom SLY token object."""
        token = Token()
        token.value = cx
        token.type = cat
        token.index = index
        token.lineno = 0
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
            
            # add custom tokens at end of word
            for ending in self.end_tokens(slot):
                yield self.get_token(ending, ending.cat, i)
                i += 1
        
    def end_tokens(self, slot):
        """Add disambiguator tokens at the end of words."""

        F = self.tf_api.F
        pos = self.slot2pos[str(slot)]
        prs = F.prs.v(slot)
        st = F.st.v(slot)

        # add state token
        if pos != 'PREP':
            if pos not in {'CARD', 'CARD1', 'ADVB'} and st == 'c':
                yield CX(slot, 'C')
            elif pos in {'NOUN','PROPN'} and prs in {'absent', 'n/a'}:
                yield CX(slot, 'A')

        # add suffix token
        if F.prs.v(slot) not in {'absent', 'n/a'}:
            yield CX(slot, 'SFX')

        # add separator token
        septoken = self.get_septoken(slot)
        if septoken:
            yield CX(slot, septoken)

    def get_septoken(self, slot):
        """Retrive separator tokens of various kinds."""
    
        if slot not in self.word_seps:
            return None

        # process various flavors of sep tokens based on context 
        P = self.getP(slot)
        if P(1,'pdp') == 'prep':
            return 'SEPPREP'
        elif P(1,'st') == 'c':
            return 'SEPGP'
        else:
            return 'SEP'
        
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
            self.error_tracker['e'] = 'reached end'

    debugfile = 'parser.out'
    
    # -- FINAL PHRASE --
    @_('np', 'adjv', 'defi', 'gp',
       'pp', 'appo', 'quant', 'num',
       'para', 'demon', 'cardc', 'advb',
       'quant_sized', 'adjv_para', 'gam_adjv',
       'gam_para',)
    def phrase(self, p):
        return p[0]

    # -- noun-based phrases --
    @_('NOUN A', 'PROPN A')
    def np(self, p):
        return p[0].slot

    @_('NOUN SFX', 'QUANT SFX', 
       'CARD1 SFX')
    def np(self, p):
        return p[0].slot

    # -- adjective-based phrases --
    @_('np ADJV', 'np ADVB', 'np ORDN') # num ADJV
    def adjv(self, p):
        return [p[1].slot, p[0], 'ADJV']

    @_('CARD1 ADJV')
    def adjv(self, p):
        return [p[1].slot, p[0].slot, 'ADJV']

    @_('np advb', 'np adjv_para')
    def adjv(self, p):
        return [p[1], p[0], 'ADJV']

    @_('ADVB PRPS', 'ADVB CARD1')
    def adjv(self, p):
        return [p[0].slot, p[1].slot, 'ADJV'] 

    @_('ADJV np', 'NEGA np')
    def adjv(self, p):
        return [p[0].slot, p[1], 'ADJV']

    @_('ADVB np', 'ADVB pp', 
       'ADVB quant', 'PRIN np',
       'ADVB num', 'ADVB defi',
       'ADVB gp')
    def adjv(self, p):
        return [p[0].slot, p[1], 'ADJV'] 

    # NB: Work-around solution for these:
    # point them at themselves, since they
    # modify their suffix; too bad the ETCBC
    # doesn't use suffix slots
    @_('ADVB SFX')
    def adjv(self, p):
        return [p[0].slot, p[0].slot, 'ADJV']

    @_('adjv ADJV')
    def adjv(self, p):
        return [p[1].slot, p[0], 'ADJV']

    @_('GAM np', 'GAM pp', 
       'GAM quant', 'GAM num', 
       'GAM defi', 'GAM gp')
    def gam_adjv(self, p):
        return [p[0].slot, p[1], 'ADJV'] 

    @_('GAM PRPS', 'GAM CARD1', 'GAM ADVB')
    def gam_adjv(self, p):
        return [p[0].slot, p[1].slot, 'ADJV'] 

    # phrase + PP modification
    @_('CARD1 adjv_pp')
    def adjv(self, p):
        return [p[1], p[0].slot, 'ADJV']

    @_('np adjv_pp')
    def adjv(self, p):
        return [p[1], p[0], 'ADJV']

    # adjectivals with construct
#    @_('ADJV C')
#    def adjconstr(self, p):
#        return p[0].slot

#    @_('adjvconstr phrase')
#    def adjv(self, p):
#        return [p[0], p[1], 'ADJV']
 
    # -- adverbial phrases -- 
    # TODO: add rule to recognize distributive 
    # constructions suhch as סביב סביב
    @_('ADJV ADVB', 'ADJV ADJV')
    def advb(self, p):
        return [p[1].slot, p[0].slot, 'ADVB']

    @_('ADVB ADVB')
    def advb(self, p):
        return [p[0].slot, p[1].slot, 'ADVB']

    # -- definite phrases --
    @_('ART ORDN', 'ART PRDE', 'ART QUANT', 
       'ART CARD1', 'ART ADJV', 
       'ART PREP', # covers some borderline preps (e.g. עֵבֶר 
      )
    def defi(self, p):
        return [p[0].slot, p[1].slot, 'DEF']

    @_('ART adjv', 'ART np', 'ART gp')
    def defi(self, p):
        return [p[0].slot, p[1], 'DEF']

    # -- genitive phrases --
    @_('NOUN C', 'PROPN C', 'ADJV C')
    def constr(self, p): 
        return p[0].slot

    @_('constr NOUN', 'constr PROPN',
       'constr QUANT', 'constr CARD1',
       'constr PRIN', 'constr ORDN',
       'constr ADVB', 'constr PRPS')
    def gp(self, p):
        return [p[1].slot, p[0], 'GP']

    @_('constr phrase')
    def gp(self, p):
        return [p[1], p[0], 'GP']

    # -- quantified phrases -- 
    @_('QUANT C', 'ORDN C')
    def qconstr(self, p):
        return p[0].slot

    @_('qconstr PROPN', 'qconstr PRDE')
    def quant(self, p):
        return [p[0], p[1].slot, 'QUANT'] 

    @_('qconstr phrase')
    def quant(self, p):
        return [p[0], p[1], 'QUANT']

    @_('QUANT np')
    def quant(self, p):
        return [p[0].slot, p[1], 'QUANT']

    @_('np QUANT')
    def quant(self, p):
        return [p[1].slot, p[0], 'QUANT']
    
    @_('QUANT PRDE')
    def quant(self, p):
        return [p[0].slot, p[1].slot, 'QUANT']

    # e.g. רב מאד
    # NB: need to add recognition for distributive cx: e.g. מעט מעט
    # NOTE: where BHSA ls==nmdi it is a "distributive noun", this feature
    # might be used to capture more distributive cases
    @_('QUANT ADVB', 'QUANT QUANT')
    def quant_sized(self, p):
        return [p[1].slot, p[0].slot, 'ADJV']

    @_('np quant_sized')
    def quant(self, p):
        return [p[1], p[0], 'QUANT']
    
    # -- demonstrative phrases --
    @_('CARD1 PRDE', 'PRIN PRDE', 'ADVB PRDE')
    def demon(self, p):
        return [p[1].slot, p[0].slot, 'DEMON']

    @_('NOUN C PRDE')
    def demon(self, p):
        return [p[2].slot, p[0].slot, 'DEMON']

    @_('np PRDE', 'defi PRDE')
    def demon(self, p):
        return [p[1].slot, p[0], 'DEMON']

    # -- appositional phrases --
    @_('np np', 'np defi', 'defi defi', 
       'appo appo', 'adjv adjv', 'gp gp',
       'appo np', 'appo defi', 'num num',
       'defi np', 'advb advb') 
    def appo(self, p):
        return [p[1], p[0], 'APPO']

    # -- prepositional phrases --
    @_('PREP np', 'PREP adjv', 'PREP defi', 'PREP gp',
       'PREP pp', 'PREP appo', 'PREP quant', 'PREP num',
       'PREP para', 'PREP demon', 'PREP cardc', 'PREP advb',
       'PREP quant_sized')
    def pp(self, p):
        return [p[0].slot, p[1], 'PP']

    @_('PREP ADVB', 'PREP PROPN', 'PREP CARD1', 
       'PREP PRDE', 'PREP QUANT', 'PREP INRG', 
       'PREP PRIN', 'PREP PRPS', 'PREP ORDN')
    def pp(self, p):
        return [p[0].slot, p[1].slot, 'PP']
    
    @_('PREP SFX')
    def pp(self, p):
        return p[0].slot

    @_('PREP PREP')
    def pp(self, p):
        return [p[0].slot, p[1].slot, 'PP']

    @_('ADJVPREP defi', 'ADJVPREP np',
       'ADJVPREP adjv', 'ADJVPREP gp',
       'ADJVPREP num')
    def adjv_pp(self, p):
        return [p[0].slot, p[1], 'PP']

    @_('ADJVPREP SFX')
    def adjv_pp(self, p):
        return p[0].slot
 
    # -- parallel phrases --
    @_('np CONJ np', 'np CONJ defi', 'np CONJ para', 
       'np CONJ adjv', 'np CONJ gp', 'np CONJ quant',
       'defi CONJ defi', 'defi CONJ para', 'defi CONJ np',
       'defi CONJ gp',
       'appo CONJ appo', 'appo CONJ para', 
       'adjv CONJ np',
       'demon CONJ demon', 'demon CONJ para', 
       'para CONJ para', 
       'num CONJ num', 'num CONJ np',
    ) 
    def para(self, p):
        conj = [p[1].slot, p[2], 'CONJ']
        return [conj, p[0], 'PARA']

    @_('np SEP np', 'np SEP para',
       'defi SEP defi', 'defi SEP para',
       'gp SEPGP gp', 'gp SEPGP para',
       'num SEP num', 'num SEP para',
       'appo SEP appo', 'appo SEP para',
       'pp SEPPREP pp', 'pp SEPPREP para',
       'para SEP para')
    def para(self, p):
        return [p[2], p[0], 'PARA']

    # Unlike, e.g. np np, a double preposition
    # pattern typically indicates parallelism
    @_('pp pp')
    def para(self, p):
        return [p[1], p[0], 'PARA']

    @_('PRPS CONJ PRPS')
    def para(self, p):
        conj = [p[1].slot, p[2].slot, 'CONJ']
        return [conj, p[0].slot, 'PARA']

    @_('PRPS CONJ np', 'PRPS CONJ quant',
       'PRPS CONJ defi', 'PRPS CONJ gp',
       'PRPS CONJ para')
    def para(self, p):
        conj = [p[1].slot, p[2], 'CONJ']
        return [conj, p[0].slot, 'PARA']

    @_('ADVB CONJ ADVB', 'ADJV CONJ ADJV')
    def adjv_para(self, p):
        conj = [p[1].slot, p[2].slot, 'CONJ']
        return [conj, p[0].slot, 'PARA']

    @_('np CONJ QUANT')
    def para(self, p):
        conj = [p[1].slot, p[2].slot, 'CONJ'] 
        return [conj, p[0], 'PARA']

    @_('pp CONJPP pp', 'pp CONJPP para',
       'gp CONJGP gp', 'gp CONJGP para',
       'adjv CONJADJV adjv',
       'quant CONJQUANT quant', 'quant CONJQUANT para')
    def para(self, p):
        conj = [p[1].slot, p[2], 'CONJ']
        return [conj, p[0], 'PARA']

    # special cases of parallelism
    @_('gam_adjv gam_adjv', 'gam_adjv gam_para')
    def gam_para(self, p):
        return [p[1], p[0], 'PARA']

    @_('gam_adjv SEP gam_adjv')
    def gam_para(self, p):
        return [p[2], p[0], 'PARA']

    # -- chained cardinal numbers --
    @_('CARD CARD', 'CARD C CARD')
    def cardc(self, p):
        return [p.CARD0.slot, p.CARD1.slot, 'CARDC']

    @_('cardc CARD')
    def cardc(self, p):
        return [p[0], p[1].slot, 'CARDC']

    @_('CARD cardc')
    def cardc(self, p):
        return [p[0].slot, p[1], 'CARDC']

    @_('CARD CONJCARD cardc')
    def cardc(self, p):
        conj = [p[1].slot, p[2], 'CONJ']
        return [conj, p[0].slot, 'CARDC']

    @_('CARD CONJCARD CARD')
    def cardc(self, p):
        conj = [p[1].slot, p[2].slot, 'CONJ']
        return [conj, p[0].slot, 'CARDC']

    # -- cardinal quantifications --
    @_('np cardc', 'gp cardc', 'appo cardc', 
       'adjv cardc')
    def num(self, p):
        return [p[1], p[0], 'NUM']
    
    @_('cardc np', 'cardc gp', 'cardc appo',
       'cardc adjv', 'cardc defi')
    def num(self, p):
        return [p[0], p[1], 'NUM']

    @_('CARD1 np', 'CARD1 gp', 'CARD1 appo',
       'CARD1 adjv', 'CARD1 defi')
    def num(self, p):
        return [p[0].slot, p[1], 'NUM']

    @_('CARD1 QUANT', 'CARD1 PRPS')
    def num(self, p):
        return [p[0].slot, p[1].slot, 'NUM']

    @_('np CARD1', 'defi CARD1', 'adjv CARD1',
       'appo CARD1')
    def num(self, p):
        return [p[1].slot, p[0], 'NUM']

def parse_phrases(paths, API):
    """Apply the parser and return metrics on unmatches."""

    # load pre-processed parts of speech values
    with open(paths['slot2pos'], 'r') as infile:
        slot2pos = json.load(infile)

    # load samples
    with open(paths['samples'], 'r') as infile:
        samples = json.load(infile)

    # load disambiguators
    with open(paths['paraseps'], 'r') as infile:
        word_seps = set(json.load(infile))

    # initialize lexer/parser
    lexer = BhsaLexer(
        API, 
        slot2pos,
        word_seps,
    )

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
        if parsing is not None and error_tracker['e'] is None:
            parsed[ph_node] = parsing
        else:
            toks = [(t.type, t.value.slot) for t in tokens]
            errors[ph_node] = (error_tracker['e'], str(toks))

            # reset error tracker
            error_tracker['e'] = None 

    # export
    with open(paths['parsed'], 'w') as outfile:
        json.dump(parsed, outfile, indent=2)

    with open(paths['notparsed'], 'w') as outfile:
        json.dump(errors, outfile, indent=2)
