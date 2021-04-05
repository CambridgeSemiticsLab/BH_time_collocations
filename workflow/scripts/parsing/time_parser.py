import json
import html
from tf.fabric import Fabric
from sly import Parser as SlyParser
from sly.lex import Token as SlyToken
from tools import nav_tree as nt
from tools.load_parse import ParseLoader

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
    
    def tokenize(self, parsedphrase, start=True):
        """Follow path to right-most item (head) and yield tokens.
        
        tokenize must adjudicate what gets yielded as a token
        and what does not. At times, we need to direct tokenize down
        a path that is not on the head-path (e.g. parallel relas).
        Sometimes whether a rela becomes a token is conditional 
        on a number of factors tuned to produce good tokens.
        """
        
        # ignore these relations
        ignore = {None, 'DEF', 'ADJV', 'GP', 'APPO'}
        def tag_rela(rela, name):
            return rela == name and rela not in ignore

        # retrieve the head path and begin walking down it
        if len(parsedphrase) > 1:
            head_phrases = list(nt.get_head_path(parsedphrase))
            if parsedphrase[-1] != 'PP' and start:
                yield self.get_sly_token('Ø', 0)
        else:
            head_phrases = [[None, parsedphrase[0], None]]
            if start:
                yield self.get_sly_token('Ø', 0)

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
                    yield from self.tokenize(subphrase, start=False) # recursively tokenize them
                
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
        'Ø',
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
    @_('atelic_ext', 'simul', 'in_dur',
       'anterior', 'anterior_dur', 'posts',
       'posterior')
    def category(self, p): 
        return p[0] 

    # -- atelic extent
    @_('Ø duration')
    def atelic_ext(self, p):
        p[1].update({
            'function': 'atelic_ext',
        })
        return p[1]

    # -- simultaneous --
    @_('B time', 'Ø time')
    def simul(self, p):
        return {
            'function': 'simultaneous',
        } 
    
    # -- telic extent / distances --
    @_('B duration')
    def in_dur(self, p):
        p[1].update({
            'function':'telic_ext, dist_fut, dist_past',
        })
        return p[1]

    # -- anterior --
    @_('L PNH/ duration', 'L PNH/ time')
    def anterior(self, p):
        p[2].update({
            'function': 'anterior',
        })
        return p[2]

    # -- anterior durative --
    @_('<D duration', '<D time')
    def anterior_dur(self, p):
        p[1].update({
            'function': 'ant_dur',
        })
        return p[1]
    
    # -- posteriors --
    @_('MN time', 'MN duration')
    def posts(self, p):
        p[1].update({
            'function': 'post, post_dur',
        })
        return p[1]
 
    # -- posterior durative --
    @_('>XR/ duration', '>XR/ time')
    def posterior(self, p):
        p[1].update({
            'function': 'posterior'
        })
        return p[1]

    # -- reanalyze "end of time" as point
    @_('QY/ duration', 'QY/ time')
    def time(self, p):
        p[1].update({
            'time': 'singular'
        })
        return p[1]
 
    # -- durations --
    @_('TIMES')
    def duration(self, p):
        return {
            'time': 'durative'
        } 

    @_('NUM duration', 'NUM time')
    def duration(self, p):
        p[1].update({
            'time': 'durative',
            'mensural': True,
        })
        return p[1]

    @_('duration duration')
    def duration(self, p):
        p[1].update(p[0])
        return p[1]

    # -- time --
    @_('TIME')
    def time(self, p):
        return {
            'time': '?',
        }
  
def parse_times(paths, API):
    """Apply the time parser."""

    # load phrase parsings
    phrases = ParseLoader(paths['ph_parses']).load()

    # initialize tokenizer
    tokenizer = TimeTokenizer(API)

    # initialize parser
    # error_tracker allows us to save
    # an offending value and keep it associated with 
    # the respective phrase node number;
    # gets reset during each iteration of the loop
    error_tracker = {'e': None}
    parser = TimeParser(error_tracker)

    # run the parser
    parsed = {}
    errors = {}
    for ph_node, parsing in phrases.items():

        # skip non-time phrases
        phrase = API.L.u(ph_node, 'phrase')[0]
        if API.F.function.v(phrase) != 'Time':
            continue

        tokens = list(tokenizer.tokenize(parsing))
        parsing = parser.parse(t for t in tokens)
        if parsing is not None and error_tracker['e'] is None:
            parsed[ph_node] = parsing
        else:
            toks = html.escape(str([t.type for t in tokens]))
            errors[ph_node] = (error_tracker['e'], toks)

            # reset error tracker
            error_tracker['e'] = None 

    # export
    with open(paths['parsed'], 'w') as outfile:
        json.dump(parsed, outfile, indent=2, ensure_ascii=False)

    with open(paths['notparsed'], 'w') as outfile:
        json.dump(errors, outfile, indent=2, ensure_ascii=False)
