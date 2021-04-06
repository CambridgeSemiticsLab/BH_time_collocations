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

    def __init__(self, paths, tf_api):
        self.F, self.L = tf_api.F, tf_api.L
        self.tokens = set() # to feed to the parser later
        self.slot2pos = ParseLoader(paths['slot2pos']).load()
        with open(paths['lexmap'], 'r') as infile:
            self.lexmap = json.load(infile)
    
    def sly_token(self, tag, index):
        """Get SLY token object with customized data."""
        token = SlyToken()
        self.tokens.add(tag) # track unique tokens
        token.value = ''
        token.type = tag
        token.index = index
        token.lineno = 1
        return token

    def get_head(self, phrase):
        """Flexible head grabbing, allowing for singular phrases."""
        if type(phrase) == int:
            return phrase
        elif len(phrase) == 1:
            return phrase[0]
        else:
            return nt.get_head(phrase)
    
    def tokenize(self, parsedphrase, start=True):
        """Follow path to right-most item (head) and yield tokens.
        
        tokenize must adjudicate what gets yielded as a token
        and what does not. At times, we need to direct tokenize down
        a path that is not on the head-path (e.g. parallel relas).
        Sometimes whether a rela becomes a token is conditional 
        on a number of factors tuned to produce good tokens.
        """
        
        # ignore these relations
        ignore = {None, 'ADJV'}
        def tag_rela(rela, name):
            return rela == name and rela not in ignore

        # retrieve the head path and begin walking down it
        if len(parsedphrase) == 3:
            head_phrases = list(nt.get_head_path(parsedphrase))
        else:
            head_phrases = [[None, parsedphrase[0], None]]

        for i, phrase in enumerate(head_phrases):
            
            src, tgt, rela = phrase

            # yield null token if no starting preposition found
            if (i==0) and start:
                null = self.null_token(phrase, i)
                if null:
                    yield self.sly_token('Ø', 0)

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
            
            # numbered relations
            elif tag_rela(rela, 'NUM'):            
                yield self.num_token(phrase, i)

            # quantified relations
            elif tag_rela(rela, 'QUANT'):
                token = self.quant_token(phrase, i)
                if token:
                    yield token

            # definite relations
            elif tag_rela(rela, 'DEF'):
                yield self.def_token(phrase, i)

            elif tag_rela(rela, 'DEMON'):
                yield self.demon_token(phrase, i)

            # genitival relations
            elif tag_rela(rela, 'GP'):
                token = self.gen_token(phrase, i)
                if token:
                    yield token

            # drip-bucket tokenizer
            elif rela not in ignore:
                yield self.sly_token(rela, i)
                
            # give TIME token once reaching the head
            if i+1 == len(head_phrases):
                yield from self.time_token(phrase, i)
                        
    def time_token(self, phrase, i):
        """Parse times into tokens."""
        time = phrase[1]

        # yield items attached to the time
        if self.F.prs.v(time) not in {'absent', 'n/a'}:
            yield self.sly_token('SFX', i)
        if self.F.nu.v(time) == 'du':
            yield self.sly_token('NUM', i)

        # yield the time itself
        lex = self.F.lex.v(time)
        if self.slot2pos[time] == 'ADVB':
            token = self.lexmap.get(lex, 'ADVB')
        elif self.slot2pos[time] == 'PREP':
            token = self.F.lex.v(time)
        elif self.F.nu.v(time) == 'pl':
            token = 'TIMES'
        else:
            token = 'TIME'
        yield self.sly_token(token, i)
    
    def prep_token(self, phrase, i):
        """Parse prepositions into singular tokens."""
        token = self.F.lex.v(phrase[0])
        token = self.lexmap.get(token, token)
        return self.sly_token(token, i)
    
    def null_token(self, phrase, i):
        """Parse whether phrase begins without preposition."""
        
        # process single-word phrases
        if phrase[-1] is None:
            if self.slot2pos[phrase[1]] != 'PREP':
                return self.sly_token('Ø', i)
            else:
                return None # prevent running rest

        # process multi-word phrases;
        # a bit complicated process; we want to look for a 
        # first word that is not a preposition, unless the
        # first word is an adverb, then we ignore it and check
        # the next first word
        slots = sorted(
            s for s in nt.get_slots(phrase)
                if self.F.pdp.v(s) != 'advb'
        )
        is_null = (
            (not slots) # only advb
            or (self.slot2pos[slots[0]] != 'PREP')
        )
        if is_null:
            return self.sly_token('Ø', i)

    def def_token(self, phrase, i):
        """Definite tokens."""
        return self.sly_token('THE', i)

    def appo_token(self, phrase, i):
        """Parse appositional tokens."""
        src, tgt, rela = phrase
        
        # get the head item of the appositional phrase
        appo_head = self.get_head(src)
        
        # process appositional demonstratives
        if self.F.pdp.v(appo_head) == 'prde':
            appo_lex = self.F.lex.v(appo_head)
            token = self.lexmap[appo_lex]
            return self.sly_token(token, i)
            
    def demon_token(self, phrase, i):
        """Parse demonstratives."""
        src, tgt, rela = phrase
        demon = self.get_head(src)
        demonlex = self.F.lex.v(demon)
        token =self.lexmap[demonlex]
        return self.sly_token(token, i)

    def num_token(self, phrase, i):
        """Parse number tokens."""
        number = phrase[0]
        is_one = (
            type(number) == int
            and self.F.lex.v(number) == '>XD/'
            and type(phrase[1]) == int
        )
        if is_one:
            token = 'NUM_ONE'
        else:
            token = 'NUM'
        return self.sly_token(token, i)

    def gen_token(self, phrase, i):
        """Parse genitive phrase tokens."""
        src, tgt, rel = phrase
        gen_head = self.get_head(src)
        ph_head = self.get_head(tgt)

        # identify genitive durations
        # e.g. חדשׁ ימים "month of days"
        dur_lexs = {'JWM/', 'DWR/', 'NYX/'}
        gen_dur = (
            self.F.nu.v(ph_head) == 'sg'
            and self.F.lex.v(gen_head) in dur_lexs
            and self.F.nu.v(gen_head) == 'pl'
        )
        if gen_dur:
            return self.sly_token('GENDUR', i)

    def quant_token(self, phrase, i):
        """Return quantifier tokens."""
        src, tgt, rel = phrase
        quant = self.get_head(src)
        quantlex = self.F.lex.v(quant)
        token = self.lexmap.get(quantlex)
        if token:
            return self.sly_token(token, i) 

class TimeParser(SlyParser):

    """Parse semantic tokens with a YACC grammar."""

    # initialize standard methods / attributes
    def __init__(self, error_tracker):
        super().__init__()
        self.error_tracker = error_tracker

    tokens = {
        'TIME',
        'TIMES',
        'Ø',
        '<D',
        '>XR/',
        'B',
        'L',
        'MN',
        'K',
        'PNH/',
        'THE', 'THIS', 'THAT',
        'SFX',
        'NUM',
        'NUM_ONE',
        'BEGINNING',
        'MIDDLE',
        'END',
        'MANY',
        'ALL',
        'GENDUR',
        'NOW',
        'TOMORROW',
        'YESTERDAY',
        'DURATION',
        'THUS',
        'THEN',
    }

    def error(self, token):
        """Keep track of errors."""
        try:
            self.error_tracker['e'] = token.value.slot
        except:
            self.error_tracker['e'] = 'reached end'

    debugfile = '../results/data_metrics/timeparse.debug'
       
    # -- FINAL MATCHES --
    @_('atelic_ext', 'simul', 'in_dur',
       'anterior', 'anterior_dur', 'posts',
       'posterior', 'atelic_simul', 'antdur_simul')
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
    @_('B time', 'Ø time',
       'K time')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
        })
        return p[1] 
    
    # either atelic ext or simul
    @_('Ø dur_sing')
    def atelic_simul(self, p):
        p[1].update({
            'function': 'atelic_ext, simultaneous',
        })
        return p[1]
    
    # -- telic extent / distances --
    @_('B duration', 'B dur_sing')
    def in_dur(self, p):
        p[1].update({
            'function':'telic_ext, dist_fut, dist_past',
        })
        return p[1]

    # -- anterior --
    @_('L PNH/ duration', 'L PNH/ time', 'L PNH/ dur_sing')
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
    
    # -- anterior dur / simultaneous
    @_('L duration', 'L time')
    def antdur_simul(self, p):
        p[1].update({
            'function': 'ant_dur, simultaneous',
        })
        return p[1]

    # -- posteriors --
    @_('>XR/ duration', '>XR/ time', '>XR/ dur_sing')
    def posterior(self, p):
        p[1].update({
            'function': 'posterior',
        })
        return p[1]

    @_('SFX >XR/')
    def posterior(self, p):
        return {
            'function': 'posterior',
            'reference': 'deictic',
            'ref_type': 'personal',
        }

    # stand-alone >XR, will only match if 
    # followed by nothing else; causes shift/reduce
    # conflict, but this is tolerable for this proj.
    @_('>XR/')
    def posterior(self, p):
        return {
            'function': 'posterior',
        }

    # -- posteriors (durative?) --
    @_('MN time', 'MN duration', '<D dur_sing')
    def posts(self, p):
        p[1].update({
            'function': 'post, post_dur',
        })
        return p[1]

    # -- reanalyze "end of time" as point
    @_('BEGINNING duration', 'BEGINNING time')
    def time(self, p):
        p[1].update({
            'time': 'singular',
            'time_loc': 'beginning',
        })
        return p[1]

    @_('MIDDLE duration', 'MIDDLE time')
    def time(self, p):
        p[1].update({
            'time': 'singular',
            'time_loc': 'middle',
        })
        return p[1]

    @_('END duration', 'END time')
    def time(self, p):
        p[1].update({
            'time': 'singular',
            'time_loc': 'end',
        })
        return p[1]
 
    # -- durations --
    @_('TIMES')
    def duration(self, p):
        return {
            'time': 'durative',
        } 

    @_('DURATION')
    def duration(self, p):
        return {
            'time': 'durative',
            'advb': True,
        }

    @_('NUM duration', 'NUM time',
       'NUM_ONE duration')
    def duration(self, p):
        p[1].update({
            'time': 'durative',
            'mensural': True,
        })
        return p[1]

    @_('ALL duration', 'ALL time',
       'MANY duration', 'MANY time')
    def duration(self, p):
        p[1].update({
            'time': 'durative',
        })
        return p[1]

    @_('duration duration')
    def duration(self, p):
        p[1].update(p[0])
        return p[1]

    @_('GENDUR time')
    def duration(self, p):
        p[1].update({
            'time': 'durative',
        })
        return p[1]

    # add reference data
    @_('SFX duration')
    def duration(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'personal',
        })
        return p[1]

    @_('THE duration')
    def duration(self, p):
        p[1].update({
            'reference': 'anaphora',
            'ref_type': 'exophoric',
        })
        return p[1]

    @_('THIS duration')
    def duration(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'near'
        })
        return p[1]

    @_('THAT duration')
    def duration(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'far'
        })
        return p[1]

    @_('NUM_ONE time')
    def dur_sing(self, p):
        p[1].update({
            'time': 'durative, singular', # disambig needed  
            'mensural': True,
        })
        return p[1]

    # -- time --
    @_('TIME')
    def time(self, p):
        return {
            'time': '?',
        }
   
    @_('SFX time')
    def time(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'personal',
        })
        return p[1]

    @_('THE time')
    def time(self, p):
        p[1].update({
            'reference': 'anaphora',
            'ref_type': 'exophoric',
        })
        return p[1]

    @_('THIS time')
    def time(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'near'
        })
        return p[1]

    @_('THAT time')
    def time(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'far'
        })
        return p[1]

    # -- adverb time --
    @_('NOW')
    def time(self, p):
        return {
            'time': 'singular',
            'reference': 'deictic',
            'tense': 'present',
        }

    @_('TOMORROW')
    def time(self, p):
        return {
            'time': 'singular',
            'reference': 'deictic',
            'tense': 'future',
        }

    @_('YESTERDAY')
    def time(self, p):
        return {
            'time': 'singular',
            'reference': 'deictic',
            'tense': 'past',
        }

    @_('THUS')
    def time(self, p):
        return {
            'time': 'singular',
            'reference': 'deictic',
            'ref_type': 'textual',
        }

    @_('THEN')
    def time(self, p):
        return {
            'time': 'singular',
            'reference': 'deictic',
            'ref_dist': 'far',
        }
 

def parse_times(paths, API):
    """Apply the time parser."""

    # load phrase parsings
    phrases = ParseLoader(paths['ph_parses']).load()

    # initialize tokenizer
    tokenizer = TimeTokenizer(paths, API)

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
