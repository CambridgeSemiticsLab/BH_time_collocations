import json
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

    def __init__(self, paths, tf_api, debug=False):
        self.F, self.L = tf_api.F, tf_api.L
        self.tokens = set() # to feed to the parser later
        self.slot2pos = ParseLoader(paths['slot2pos']).load()
        self.debug = debug
        with open(paths['lexmap'], 'r') as infile:
            self.lexmap = json.load(infile)
    
    def sly_token(self, tag, index, value=''):
        """Get SLY token object with customized data."""
        token = SlyToken()
        self.tokens.add(tag) # track unique tokens
        token.value = value
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

    def get_head_phrases(self, phrase):
        """Retrieve all phrases with head in the path."""
        # retrieve the head path and begin walking down it
        if len(phrase) == 3:
            head_phrases = list(nt.get_head_path(phrase))
        else:
            head_phrases = [[None, phrase[0], None]]
        return head_phrases

    def tokenize(self, parsedphrase, start=True):
        """Follow path to right-most item (head) and yield tokens.
        
        tokenize must adjudicate what gets yielded as a token
        and what does not. At times, we need to direct tokenize down
        a path that is not on the head-path (e.g. parallel relas).
        Sometimes whether a rela becomes a token is conditional 
        on a number of factors tuned to produce good tokens.
        """
        
        # ignore these relations
        ignore = {None, 'ADJV', 'ADVB', 'SPEC'}
        def tag_rela(rela, name):
            return rela == name and rela not in ignore
   
        # process all internal relations
        head_phrases = self.get_head_phrases(parsedphrase)
        paras = []
        
        i = 0
        for phrase in head_phrases:
            
            src, tgt, rela = phrase

            # yield null token if no starting preposition found
            if (i==0) and start:
                token = self.null_token(phrase, i)
                if token:
                    yield token
                    i += 1

            # tokenize prepositions
            if tag_rela(rela, 'PP'):
                yield self.prep_token(phrase, i)
                
            # handle parallels
            elif tag_rela(rela, 'PARA'):
                for subphrase in nt.unfold_paras(src):
                    paras.extend(
                        self.tokenize(subphrase, start=False) # recursively tokenize them
                    )                    
            
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

            # card chain relas:
            # NB: that this will only ever happen if 
            # a cardinal chain itself is profiled, without
            # a quantified element
            elif tag_rela(rela, 'CARDC'):
                yield self.sly_token('CARD', i)
    
            # drip-bucket tokenizer
            elif rela not in ignore:
                yield self.sly_token(rela, i)

            # increase index
            i += 1
                
            # give TIME token once reaching the end
            if head_phrases[-1] == phrase:
                timetoks = list(self.time_token(phrase, i))
                seenlex = set(tok.value for tok in timetoks)
                for tok in timetoks:
                    yield tok
                i += 1

                # yield all of the parallel tokens
                for sp_token in paras:
                    sp_token.index = i

                    # detect repeated lexemes (distributive construction)
                    if (sp_token.type in {'TIME'}
                    and sp_token.value in seenlex):
                        sp_token.type = 'TIME2'

                    yield sp_token
                    i += 1
                        
    def time_token(self, phrase, i):
        """Parse times into tokens."""
        src, time, rela = phrase

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
        elif self.slot2pos[time] == 'CARD' and rela == 'CARDC':
            token = 'CARD'
        elif self.F.nu.v(time) == 'pl':
            token = 'TIMES'
        else:
            token = 'TIME'
        yield self.sly_token(token, i, value=lex)
    
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
        head = nt.get_head(phrase)
        slots = sorted(
            s for s in nt.get_slots(phrase)
                if (self.F.pdp.v(s) != 'advb' or s == head)
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
        token = self.lexmap[demonlex]
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
        'TIME2',
        'TIMES',
        'Ø',
        '<D',
        '>XR/',
        'B',
        'L',
        'MN',
        'K',
        'BJN/',
        '>L',
        '<L',
        '>T',
        'VRM/',
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
        'CARD',
        'ONWARD',
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
       'posterior', 'atelic_simul', 'antdur_simul',
       'hab_simul', 'begin_to_end', 'multi_simul',
       'antpost_dur', 'posterior_dur', 'multi_postantdur',
       'multi_antpost', 'multi_antdursim', 'dur_to_end',
       'simul_posts', 'multi_posts', 'simul_to_end')
    def category(self, p): 
        return p[0] 

    # -- composed constructions
    @_('simul simul', 'in_dur in_dur')
    def hab_simul(self, p):
        data = {
            'function': 'simul_habitual, multi_simul',
            'parts': [p[0], p[1]]
        }
        if p[1].get('duplicate'):
            data['distributive'] = True
        return data

    @_('hab_simul in_dur')
    def hab_simul(self, p):
        data = {
            'function': 'simul_habitual, multi_simul',
            'parts': p[0]['parts'] + [p[1]]
        }
        return data

    @_('hab_simul hab_simul')
    def hab_simul(self, p):
        p[0]['parts'].extend(
            p[1]['parts']
        )
        return p[0]

    @_('posts anterior_dur', 'posts antdur_simul',
       'posterior_dur anterior_dur', 'antpost_dur anterior_dur')
    def begin_to_end(self, p):
        p[0]['function'] = 'posterior_dur'
        getattr(p, 'posts', {})['function'] = 'anterior_dur'
        data = {
            'function': 'begin_to_end',
            'parts': [p[0], p[1]]
        }
        return data

    @_('posts ONWARD')
    def begin_to_end(self, p):
        p[0]['function'] = 'posterior_dur'
        onward = {'time': 'durative'}
        data = {
            'function': 'begin_to_end',
            'parts': [p[0], onward]
        }
        return data

    # 'all the days of Koresh and on...'
    @_('atelic_ext anterior_dur')
    def dur_to_end(self, p):
        data = {
            'function': 'dur_to_end',
            'parts': [p[0], p[1]]
        }
        return data

    # 'on the eigth day and onward'
    @_('simul ONWARD')
    def simul_to_end(self, p):
        data = {
            'function': 'simul_to_end',
            'parts': [p[0], {'time': 'single', 'advb': True}]
        }
        return data

    @_('in_dur simul', 'simul in_dur')
    def multi_simul(self, p):
        p.in_dur['function'] = 'telic_ext'
        data = {
            'function': 'multi_simul',
            'parts': [p[0], p[1]]
        }
        return data

    @_('multi_simul simul')
    def multi_simul(self, p):
        parts = p[0]['parts'] + [p[1]]
        data = {
            'function': 'multi_simul',
            'parts': parts,
        }
        return data

    @_('posts anterior')
    def multi_postantdur(self, p):
        p[0]['function'] = 'anterior_dur'
        data = {
            'function': 'multi_postant',
            'parts': [p[0], p[1]]
        }
        return data

    @_('anterior posterior')
    def multi_antpost(self, p):
        data = {
            'function': 'multi_antpost',
            'parts': [p[0], p[1]],
        }
        return data

    @_('antdur_simul antdur_simul')
    def multi_antdursim(self, p):
        data = {
            'function': 'multi_simul, multi_anterior_dur',
            'parts': [p[0], p[1]]
        }
        return data

    @_('antdur_simul multi_antdursim')
    def multi_antdursim(self, p):
        partsb = p[1]['parts'] 
        p[1]['parts'] = [p[0]] + partsb
        return p[1]

    # 'in that day and after tomorrow'
    @_('simul posts')
    def simul_posts(self, p):
        data = {
            'function': 'simul_posts',
            'parts': [p[0], p[1]]
        }
        return data

    @_('posts posts')
    def multi_posts(self, p):
        data = {
            'function': 'multi_posts',
            'parts': [p[0], p[1]]
        }
        return data

    # -- atelic extent
    @_('Ø duration')
    def atelic_ext(self, p):
        p[1].update({
            'function': 'atelic_ext',
        })
        return p[1]

    # NB: this pattern is very significant
    # as it corroborates Haspelmath's hypothesis
    # that the zero-marked atelic extent derives
    # from the accusative / object role
    @_('>T duration')
    def atelic_ext(self, p):
        p[1].update({
            'function': 'atelic_ext',
        })
        return p[1]

    @_('atelic_ext atelic_ext')
    def atelic_ext(self, p):
        data = {
            'function': 'atelic_ext',
            'parts': [p[0], p[1]]
        }
        return data

    # -- simultaneous --
    @_('B time', 'Ø time',
       'K time', 'BJN/ duration',
       '>L time', '<L time')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
        })
        return p[1] 

    @_('B simul', 'K simul', 'L simul')
    def simul(self, p):
        return p[1]

    # Very interesting pattern here, somewhat
    # similar to antpost_dur:
    # למן־היום (II Sam 19:25)
    # except here I think this is a simultaneous
    # position after a duration: 
    # כמשלש חדשים (Gen 38:24)
    @_('K posterior_dur')
    def simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': [
                {'function': 'simultaneous'},
                p[1],
            ]
        } 
        return data

    @_('BEGINNING duration', 'BEGINNING time')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
            'time': 'singular',
            'time_loc': 'beginning',
        })
        return p[1]

    @_('MIDDLE duration', 'MIDDLE time')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
            'time': 'singular',
            'time_loc': 'middle',
        })
        return p[1]

    @_('Ø MIDDLE time')
    def simul(self, p):
        p[2].update({
            'function': 'simultaneous',
            'time': 'singular',
            'time_loc': 'middle',
        })
        return p[2]
 
    @_('END duration', 'END time')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
            'time': 'singular',
            'time_loc': 'end',
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
    @_('B duration', 'B dur_sing', 'K duration')
    def in_dur(self, p):
        p[1].update({
            'function':'telic_ext, dist_fut, dist_past',
        })
        return p[1]

    # -- anterior --
    @_('L PNH/ duration', 'L PNH/ time', 
       'L PNH/ dur_sing')
    def anterior(self, p):
        p[2].update({
            'function': 'anterior',
        })
        return p[2]

    @_('L SFX PNH/')
    def anterior(self, p):
        data = {
            'function': 'anterior',
            'reference': 'deictic',
            'ref_type': 'personal',
        }
        return data

    @_('L PNH/ posts')
    def anterior(self, p):
        data = {
            'function': 'anterior',
            'parts': [
                {'function': 'anterior'},
                p[2]
            ]
        }
        return data

    # stand-alone >XR, will only match if 
    # followed by nothing else; causes shift/reduce
    # conflict, but this is tolerable for this proj.
    @_('VRM/')
    def anterior(self, p):
        return {
            'function': 'anterior',
            'advb': True,
        }

    # -- anterior durative --
    @_('<D duration', '<D time', '<D simul',
       'L posterior', '<D dur_sing')
    def anterior_dur(self, p):
        p[1].update({
            'function': 'anterior_dur',
        })
        return p[1]

    # -- anterior dur / simultaneous
    @_('L time', 'L duration')
    def antdur_simul(self, p):
        p[1].update({
            'function': 'anterior_dur, simultaneous',
        })
        return p[1]

    @_('L dur_sing')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
            'time': 'singular',
        })
        return p[1]

    @_('<D antdur_simul')
    def anterior_dur(self, p):
        p[1].update({
            'function': 'anterior_dur'
        })
        return p[1]

    @_('L posts', '<D posts', 'L posterior_dur')
    def antpost_dur(self, p):
        p[1].update({
            'function': 'posterior_dur',   
        })
        data = {
            'function': 'antpost_dur',
            'parts': [
                {'function': 'anterior_dur'},
                p[1]
            ]
        }
        return data

    @_('<D posterior')
    def antpost_dur(self, p):
        data = {
            'function': 'antpost_dur',
            'parts': [
                {'function': 'anterior_dur'},
                p[1],
            ]
        }
        return data

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
            'advb': True,
        }

    # -- posteriors (durative?) --
    @_('MN time', 'MN dur_sing')
    def posts(self, p):
        p[1].update({
            'function': 'posterior, posterior_dur',
        })
        return p[1]

    @_('MN duration')
    def posterior_dur(self, p):
        p[1].update({
            'function': 'posterior_dur',    
        })
        return p[1]

    @_('MN simul')
    def posterior(self, p):
        p[1]['function'] = 'posterior'
        return p[1]

    @_('MN posterior')
    def posterior(self, p):
        p[1]['function'] = 'posterior'
        return p[1]

    # NB how the opposing duration leads
    # to a re-evaluation of this type!
    @_('MN anterior_dur')
    def posterior_dur(self, p):
        p[1]['function'] = 'posterior_dur'
        return p[1]

    @_('MN antdur_simul')
    def posterior_dur(self, p):
        p[1]['function'] = 'posterior_dur'
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

    @_('CARD')
    def duration(self, p):
        return {
            'time': 'durative',
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
        data = {
            'time': 'durative',
            'parts': [p[0], p[1]]   
        }
        return data

    @_('duration time', 'time duration')
    def duration(self, p):
        p.time['time'] = 'durative' # reanalyze as duration
        data = {
            'time': 'durative',
            'parts': [p[0], p[1]] 
        }
        return data 

    @_('time time')
    def duration(self, p):
        data = {
            'time': 'durative',
            'parts': [p[0], p[1]]
        }
        if p[1].get('duplicate'):
            data['distributive'] = True
        return data

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

    @_('TIME2')
    def time(self, p):
        return {
            'time': '?',
            'duplicate': True,
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
        if API.F.function.v(ph_node) != 'Time':
            continue

        debug = False
        if ph_node == 1144822:
            tokenizer.debug = True

        tokens = list(tokenizer.tokenize(parsing))
    
        parsing = parser.parse(t for t in tokens)
        if parsing is not None and error_tracker['e'] is None:
            parsed[ph_node] = parsing
        else:
            toks = str([t.type for t in tokens])
            errors[ph_node] = (error_tracker['e'], toks)

            # reset error tracker
            error_tracker['e'] = None 

    # export
    with open(paths['parsed'], 'w') as outfile:
        json.dump(parsed, outfile, indent=2, ensure_ascii=False)

    with open(paths['notparsed'], 'w') as outfile:
        json.dump(errors, outfile, indent=2, ensure_ascii=False)
