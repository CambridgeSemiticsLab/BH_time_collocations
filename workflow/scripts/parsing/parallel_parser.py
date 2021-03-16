import json
from sly import Lexer, Parser
from sly.lex import Token

class phraseLexer(Lexer):
    """Wrap parsed phrases into Token objects for further parsing."""

    # initialize standard methods / attributes
    def __init__(self):
        super().__init__()
 
    tokens = {
        ADJV,
        ADVB,
        APPO,
        CARDC,
        DEF,
        DEMON,
        GP,
        NP,
        NUM,
        PP,
        QUANT 
    }
    
    def tokenize(self, parselists):
        """Convert a bunch of parsings into tokens."""
        for i, phrase in enumerate(parselists):
            token = Token()
            token.value = phrase
            token.type = phrase[-1]
            token.index = i
            token.lineno = i
            yield token
       
class ParallelParser(Parser):

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

    debugfile = 'parser.out'
