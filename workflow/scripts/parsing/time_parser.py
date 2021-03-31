import json
from tf.fabric import Fabric
from sly import Parser as SlyParser
from sly.lex import Token as SlyToken
from scripts.tools import nav_tree as nt

class TimeTokenizer:
    """Generate selected tokens from a time phrase parsing.
    
    The tokens will be subsequently parsed to generate
    semantic function labels for the phrases.
    """

    def __init__(self, tf_api):
        self.F, self.L = tf_api.F, tf_api.L
        self.tokens = set() # to feed to the parser later
    
    def get_sly_token(self, tag, index):
        """Get SLY token object with customized data."""
        token = SlyToken()
        self.tokens.add(tag) # track unique tokens
        token.value = ''
        token.type = tag
        token.index = index
        token.lineno = 1
        return token
    
    def tokenize(self, parsedphrase):
        """Follow path to right-most item (head) and yield tokens.
        
        tokenize must adjudicate what gets yielded as a token
        and what does not. At times, we need to direct tokenize down
        a path that is not on the head-path (e.g. parallel relas).
        Sometimes whether a rela becomes a token is conditional 
        on a number of factors tuned to produce good tokens.
        """
        
        # ignore these relations
        ignore = {'DEF', 'ADJV', 'GP', 'APPO'}
        
        def tag_rela(rela, name):
            return rela == name and rela not in ignore
        
        # retrieve the head path and begin walking down it
        head_phrases = list(nt.get_head_path(parsedphrase))
        for i, phrase in enumerate(head_phrases):
            
            src, tgt, rela = phrase
            
            # tokenize prepositions
            if tag_rela(rela, 'PP'):
                yield self.prep_token(phrase, i)
                
            # handle parallels
            elif tag_rela(rela, 'PARA'):
                # we use unfold_paras in case there are 
                # additional recursively embedded parallel
                # phrases below; unfold_paras will stop as 
                # soon as it encounters a phrase with a rela
                # that != PARA or CONJ
                for subphrase in nt.unfold_paras(src):
                    yield from self.tokenize(subphrase) # recursively tokenize them
                
            # tokenize appositional relations if relevant
            elif tag_rela(rela, 'APPO'):
                token = self.appo_token(phrase, i)
                if token:
                    yield token
            
            # drip-bucket tokenizer
            elif rela not in ignore:
                yield self.get_sly_token(rela, i)
                
            # give TIME token once reaching the head
            if i+1 == len(head_phrases):
                yield self.time_token(phrase, i)
                        
    def time_token(self, phrase, i):
        """Parse times into tokens."""
        time = phrase[1]
        if self.F.nu.v(time) == 'pl':
            token = 'TIMES'
        else:
            token = 'TIME'
        return self.get_sly_token(token, i)
    
    def prep_token(self, phrase, i):
        """Parse prepositions into singular tokens."""
        token = self.F.lex.v(phrase[0])
        return self.get_sly_token(token, i)
    
    def appo_token(self, phrase, i):
        """Parse appositional tokens."""
        src, tgt, rela = phrase
        
        # get the head item of the appositional phrase
        if type(src) == int:
            appo_head = src
        else:
            appo_head = nt.get_head(phrase[0])
        
        # process appositional demonstratives
        if self.F.pdp.v(appo_head) == 'prde':
            appo_lex = self.F.lex.v(appo_head)
            dist = self.demon_dist(appo_lex).upper()
            token = f'DEM_{dist}'
            return self.get_sly_token(token, i)
            
    def demon_dist(self, lex):
        """Get the distance of a demonstrative."""
        demon_map = { 
                'Z>T': 'near',
                'HJ>': 'far',
                'HMH': 'far',
                '>LH': 'near',
                'HM': 'far',
                'HW>': 'far',
                'ZH': 'near'
        }
        return demon_map[lex]

class TimeParser(SlyParser):

    """Parse semantic tokens with a YACC grammar."""

    # initialize standard methods / attributes
    def __init__(self, error_tracker):
        super().__init__()
        self.error_tracker = error_tracker

    tokens = {
        '<D',
        '>XR/',
        'B',
        'L',
        'MN',
        'NUM',
        'PNH/',
        'QY/',
        'TIME',
        'TIMES'
    }

    def error(self, token):
        """Keep track of errors."""
        try:
            self.error_tracker['e'] = token.value.slot
        except:
            self.error_tracker['e'] = 'reached end'

    #debugfile = 'parser.out'
       
    # -- FINAL MATCHES --
    @_('simul', 'ante', 'post', 'ant_dur',
       'post_dur', 'in_dur', 'atelic_ext')
    def category(self, p): 
        return p[0] 

    @_('B TIME')
    def simul(self, p):
        return 'simul'
    
    @_('before TIME', 'before duration')
    def ante(self, p):
        return 'ante'
    
    @_('>XR/ TIME', '>XR/ duration')
    def post(self, p):
        return 'post'
    
    @_('<D TIME', '<D duration', '<D time')
    def ant_dur(self, p):
        return 'ant_dur'
    
    @_('MN TIME', 'MN time', 'MN duration')
    def post_dur(self, p):
        return 'post_dur'
    
    @_('B duration')
    def in_dur(self, p):
        return 'telic_ext|dist_fut'
    
    @_('TIME', 'duration')
    def atelic_ext(self, p):
        return 'atelic_ext'
    
    # -- TEMPORARY MATCHES --
    @_('TIMES', 'duration duration',
       'NUM TIME', 'NUM time', 'NUM duration')
    def duration(self, p):
        return ''
    
    @_('L PNH/')
    def before(self, p):
        return ''
    
    @_('QY/ duration', 'QY/ TIME')
    def time(self, p):
        return ''

def parse_times(samp_path, parsepath, noparsepath, datalocs, API=None):
    """Apply the time parser."""

    # initialize TF
    if API is None:
        TF = Fabric(locations=datalocs['bhsadata'])
        API = TF.load('lex pdp num')

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
        if parsing is not None and error_tracker['e'] is None:
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
