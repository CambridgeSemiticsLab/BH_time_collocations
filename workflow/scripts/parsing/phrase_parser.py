from sly import Lexer, Parser
from sly.lex import Token

class BhsaLexer(Lexer):
    """Wrap a BHSA node in SLY Token object.
        
    This enables us to use our custom tokens rather
    than raw string input.
    """

    # unique parts of speech values harvested from 
    # the parts of speech parser and pasted here
    # see results/parsing/uniquepos.txt
    tokens = {
        ADJV,
        ADVB,
        ART,
        CARD,
        CONJ,
        INRG,
        INTJ,
        NEGA,
        NOUN,
        ORDN,
        PRDE,
        PREP,
        PRIN,
        PROPN,
        PRPS,
        QUANT,
        VERB,
    }

    def tokenize(self, slots, slot2pos):
        """Write over the SLY tokenize method to yield custom data.

        Args:
            slots: iterable of integers which correspond with
                slots / word nodes in BHSA
            slot2pos: a dict that maps slots to their parts of speech
        
        Yields:
            generator of SLY Tokens
        """
        for i, slot in enumerate(slots):
            token = Token()
            token.value = slot
            token.type = slot2pos[slot]
            token.index = i
            token.lineno = i
            yield token

    def tokenize_morphs(self, slot):
        """Construct a token for morphological features."""
         

class PhraseParser(Parser):

    """A CFG parser for Biblical Hebrew phrases."""
    
    tokens = BhsaLexer.tokens

    @_('np', 'df', 'pp')
    def phrase(self, p):
        return p[0]

    @('ART NOUN', 'ART ADJV', 'ART ADVB')
    def df(self, p):
        return [p.ART, [p[1], 'NP']]

    @_('NOUN')
    def np(self, p):
        return [p.NOUN, 'NP']

    @_('PREP NOUN', 'PREP ADJV',
       'PREP ADVB', 'PREP df', 'PREP pp')
    def pp(self, p):
        return [p[0], p[1], 'PP']
