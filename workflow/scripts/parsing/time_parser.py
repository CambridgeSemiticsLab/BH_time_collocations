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
        self.i = 0 # token indexer

        with open(paths['lexmap'], 'r') as infile:
            self.lexmap = json.load(infile)
    
    def sly_token(self, tag, node=0, **values):
        """Get SLY token object with customized data."""
        token = SlyToken()
        self.tokens.add(tag) # track unique tokens
        token.value = dict(values)
        token.value['node'] = node
        token.value['phatom'] = self.L.u(node, 'phrase_atom')[0]
        token.type = tag
        token.index = self.i
        token.lineno = 1
        self.i += 1
        return token

    def get_head(self, phrase):
        """Flexible head grabbing, allowing for singular phrases."""
        if type(phrase) == int:
            return phrase
        elif len(phrase) == 1:
            return phrase[0]
        else:
            return nt.get_head(phrase)

    def get_slots(self, phrase):
        """Get slots without errors."""
        if type(phrase) == list and None in phrase:
            for sp in phrase:
                if type(sp) == int:
                    yield sp
        else:
            yield from nt.get_slots(phrase) 

    def get_head_phrases(self, phrase):
        """Retrieve all phrases with head in the path."""
        # retrieve the head path and begin walking down it
        if len(phrase) == 3:
            head_phrases = list(nt.get_head_path(phrase))
        else:
            head_phrases = [phrase]
        return head_phrases

    def get_rela(self, phrase):
        """Return the edge relation in a phrase."""
        if len(phrase) == 3:
            return phrase
        elif len(phrase) == 1:
            return (None, phrase[0], None)

    def tokenize(self, parsedphrase, start=True):
        """Follow path to right-most item (head) and yield tokens.
        
        tokenize must adjudicate what gets yielded as a token
        and what does not. At times, we need to direct tokenize down
        a path that is not on the head-path (e.g. parallel relas).
        Sometimes whether a rela becomes a token is conditional 
        on a number of factors tuned to produce good tokens.
        """
        
        # ignore these relations
        ignore = {None, 'ADJV', 'ADVB', 'CARDC'}
        def tag_rela(rela, name):
            return rela == name and rela not in ignore
   
        # process all internal relations
        head_phrases = self.get_head_phrases(parsedphrase)
        relas = set(self.get_rela(ph)[-1] for ph in head_phrases)
        paras = []
        
        # reset token indexer
        if start:
            self.i = 0

        for phrase in head_phrases:

            src, tgt, rela = self.get_rela(phrase) 

            # yield null token if no starting preposition found
            if (self.i==0) and start:
                token = self.null_token(phrase)
                if token:
                    yield token

            # tokenize prepositions
            if tag_rela(rela, 'PP'):
                yield self.prep_token(phrase)
                
            # handle parallels
            elif tag_rela(rela, 'PARA'):
                for subphrase in nt.unfold_paras(src):
                    paras.extend(
                        self.tokenize(subphrase, start=False) # recursively tokenize them
                    )                    
            
            # handle 'spec'
            elif tag_rela(rela, 'SPEC'):
                if type(src) == int:
                    src = [None, src, None]
                paras.extend(
                    self.tokenize(src, start=False)
                )

            # tokenize appositional relations if relevant
            elif tag_rela(rela, 'APPO'):
                token = self.appo_token(phrase)
                if token:
                    yield token
            
            # numbered relations
            elif tag_rela(rela, 'NUM'):            
                yield self.num_token(phrase)

            # quantified relations
            elif tag_rela(rela, 'QUANT'):
                token = self.quant_token(phrase)
                if token:
                    yield token

            # definite relations
            elif tag_rela(rela, 'DEF'):
                yield self.def_token(phrase)

            elif tag_rela(rela, 'DEMON'):
                yield self.demon_token(phrase)

            # genitival relations
            elif tag_rela(rela, 'GP'):
                token = self.gen_token(phrase)
                if token:
                    yield token

           # drip-bucket tokenizer
            elif rela not in ignore:
                srchead = self.get_head(src)
                yield self.sly_token(rela, node=srchead)

            # give TIME token once reaching the end
            if head_phrases[-1] == phrase:
                timetoks = list(self.time_token(phrase, relas))
                seenlex = set(tok.value.get('lex') for tok in timetoks)
                for tok in timetoks:
                    yield tok

                # yield all of the parallel tokens in sorted order
                sort_paras= lambda p: [p.value['phatom'], p.index]
                for sp_token in sorted(paras, key=sort_paras):

                    # detect repeated lexemes (distributive construction)
                    if (sp_token.type in {'TIME'}
                    and sp_token.value['lex'] in seenlex):
                        sp_token.type = 'TIME2'

                    yield sp_token
                        
    def time_token(self, phrase, time_relas):
        """Parse times into tokens."""
        src, time, rela = self.get_rela(phrase)
        F = self.F

        # yield items attached to the time
        if (prs := F.prs.v(time)) not in {'absent', 'n/a'}:
            if F.lex.v(time) == 'B' and prs == 'W':
                yield self.sly_token('IT', node=time) # בו
            else:
                yield self.sly_token('SFX', node=time)
        if F.nu.v(time) == 'du':
            yield self.sly_token('NUM', node=time)
        if F.uvf.v(time) == 'H':
            yield self.sly_token('LOCALE', node=time)

        # get data for processing time words
        calend_times = {
            'XDC=/': 'MONTH', 'CNH/': 'YEAR',
            'KSLW/': 'MONTH', '>LWL/': 'MONTH',
        }
        lex = F.lex.v(time)
        nu = F.nu.v(time)
        st = F.st.v(time)

        # process the time word itself 
        if lex in self.lexmap:
            token = self.lexmap[lex]
        elif (lex in calend_times and nu == 'sg'):
            token = calend_times[lex]
        elif lex == 'JWM/' and nu == 'sg':
            if 'NUM' in time_relas:
                token = 'DAY'
            elif (st == 'c'
                      and self.slot2pos[time+1] 
                      in {'CARD', 'CARD1'}): 
                token = 'DAY'
            elif (st == 'c'
                      and F.lex.v(time+1) == 'H'
                      and F.lex.v(time+2) == 'XDC=/'):
                token = 'DAY'
            else:
                token = 'TIME'
        elif self.slot2pos[time] == 'ADVB':
            token = 'ADVB'
        elif self.slot2pos[time] == 'PREP':
            token = F.lex.v(time)
        elif self.slot2pos[time] in {'CARD', 'CARD1'}:
            token = 'CARD'
        elif self.slot2pos[time] == 'ORDN':
            token = 'ORDN'
        elif F.nametype.v(time) == 'pers':
            token = 'PERSON'
        elif nu == 'pl':
            token = 'TIMES'
        else:
            token = 'TIME'
        yield self.sly_token(token, node=time, lex=lex)
    
    def prep_token(self, phrase):
        """Parse prepositions into singular tokens."""
        token = self.F.lex.v(phrase[0])
        token = self.lexmap.get(token, token)
        return self.sly_token(token, node=phrase[0])
    
    def null_token(self, phrase):
        """Parse whether phrase begins without preposition."""
      
        # process single-word phrases
        if len(phrase) == 1:
            if self.slot2pos[phrase[0]] != 'PREP':
                return self.sly_token('Ø', node=phrase[0])
            else:
                return None # prevent running rest

        # process multi-word phrases;
        # a bit complicated process; we want to look for a 
        # first word that is not a preposition, unless the
        # first word is an adverb, then we ignore it and check
        # the next first word
        head = nt.get_head(phrase)
        slots = sorted(nt.get_slots(phrase))
        slotsnoadvb = [
            s for s in slots
                if (
                    self.F.pdp.v(s) != 'advb'
                    or s == head
                )
        ]
        is_null = (
            (not slotsnoadvb) 
            or (self.slot2pos[slotsnoadvb[0]] != 'PREP')
        )
        if is_null:
            return self.sly_token('Ø', node=slots[0])

    def def_token(self, phrase):
        """Definite tokens."""
        return self.sly_token('THE', node=phrase[0])

    def map_demon(self, demonlex):
        demon_map = {
            "Z>T": "THIS",
            "HJ>": "THAT",
            "HMH": "THAT",
            ">LH": "THIS",
            "HM": "THAT",
            "HW>": "THAT",
            "ZH": "THIS",
        }
        return demon_map[demonlex]

    def appo_token(self, phrase):
        """Parse appositional tokens."""
        src, tgt, rela = self.get_rela(phrase)
        
        # get the head item of the appositional phrase
        appo_head = self.get_head(src)
        
        # process appositional demonstratives
        if self.F.pdp.v(appo_head) == 'prde':
            appo_lex = self.F.lex.v(appo_head)
            token = self.map_demon(appo_lex)
            return self.sly_token(token, node=appo_head)
            
    def demon_token(self, phrase):
        """Parse demonstratives."""
        src, tgt, rela = self.get_rela(phrase)
        demon = self.get_head(src)
        demonlex = self.F.lex.v(demon)
        token = self.map_demon(demonlex)
        return self.sly_token(token, node=demon)

    def num_token(self, phrase):
        """Parse number tokens."""
        number = phrase[0]
        head = self.get_head(phrase[0])
        is_one = (
            type(number) == int
            and self.F.lex.v(number) == '>XD/'
            and type(phrase[1]) == int
        )
        if is_one:
            token = 'NUM_ONE'
        else:
            token = 'NUM'
        return self.sly_token(token, node=head)

    def gen_token(self, phrase):
        """Parse genitive phrase tokens."""
        src, tgt, rel = self.get_rela(phrase)
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
            return self.sly_token('GENDUR', node=gen_head)

        # identify genitive cardinals
#        src0 = sorted(nt.get_slots(src))[0]
#        if self.slot2pos[gen_head] in {'CARD', 'CARD1'}:
#            return self.sly_token('GENCARD', node=gen_head)
#        elif self.slot2pos[src0] in {'CARD', 'CARD1'}:
#            return self.sly_token('GENCARD', node=src0)

    def quant_token(self, phrase):
        """Return quantifier tokens."""
        src, tgt, rel = self.get_rela(phrase)
        quant = self.get_head(src)
        quantlex = self.F.lex.v(quant)
        token = self.lexmap.get(quantlex)
        if token:
            return self.sly_token(token, node=quant) 

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
        'GENCARD',
        'NOW',
        'TOMORROW',
        'YESTERDAY',
        'DURATION',
        'THUS',
        'THEN',
        'CARD',
        'ONWARD',
        'MONTH',
        'YEAR',
        'DAY',
        'PERSON',
        'IT',
        'LOCALE',
        'ORDN',
        'FIRST',
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
       'simul_posts', 'multi_posts', 'simul_to_end', 
       'dur_simul', 'month_ref', 'month_simul', 'year_simul',
       'cal_simul', 'oneday_simul', 'habitual', 'day_simul',
       'ordn_simul', 'multi', 'first_simul', 'simul_ref',
       'multi_begintoend', 'posterior_dist', 'anterior_dist',
       'simul_dur', 'postdur_dist', 'antdur_dur')
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

    @_('hab_simul in_dur', 'simul hab_simul')
    def hab_simul(self, p):
        data = {
            'function': 'simul_habitual, multi_simul',
            'parts': list(p),
        }
        return data

    @_('hab_simul hab_simul')
    def hab_simul(self, p):
        data = {
            'function': 'simul_habitual, multi_simul',
            'parts': list(p),
        }
        return p[0]

    @_('Ø CARD antdur_simul')
    def hab_simul(self, p):
        p[2]['function'] = 'simultaneous'
        data = {
            'function': 'simul_habitual, multi_simul',
            'parts': [
                {'time': 'single', 'count': 'calendrical'},
                p[2]
            ]
        }
        return data

    @_('Ø year year_simul', 'Ø day day_simul')
    def habitual(self, p):
        data = {
            'function': 'habitual',
            'parts': [p[1], p[2]],
        }
        return data

    @_('posts anterior_dur', 'posts antdur_simul',
       'posterior_dur anterior_dur', 'antpost_dur anterior_dur')
    def begin_to_end(self, p):
        p[0]['function'] = 'posterior_dur'
        getattr(p, 'posts', {})['function'] = 'posterior_dur'
        getattr(p, 'antdur_simul', {})['function'] = 'anterior_dur'
        data = {
            'function': 'begin_to_end',
            'parts': [p[0], p[1]]
        }
        return data

    @_('posts ONWARD', 'posts LOCALE ONWARD', 
       'posterior_dur LOCALE duration')
    def begin_to_end(self, p):
        p[0]['function'] = 'posterior_dur'
        onward = {'time': 'durative'}
        data = {
            'function': 'begin_to_end',
            'parts': [p[0], onward]
        }
        return data

    @_('MN month month_ref')
    def begin_to_end(self, p):
        data = {
            'function': 'begin_to_end',
            'parts': [
                {'time': p[1]['time'], 'function': 'posterior_dur'},
                p[1]
            ]
        }
        return data

    @_('begin_to_end begin_to_end')
    def multi_begintoend(self, p):
        data = {
            'function': 'multi_begin_to_end',
            'parts': list(p),
        }
        return data

    @_('atelic_ext anterior_dur')
    def dur_to_end(self, p):
        data = {
            'function': 'dur_to_end',
            'parts': [p[0], p[1]]
        }
        return data

    # 'on the eigth day and onward'
    @_('simul ONWARD', 'in_dur anterior_dur', 'simul anterior_dur')
    def simul_to_end(self, p):
        data = {
            'function': 'simul_to_end',
            'parts': [p[0], {'time': 'single', 'advb': True}]
        }
        return data

    @_('in_dur simul', 'simul in_dur',
       'day_simul day_simul', 'day_simul simul',
       'year_simul simul', 'ordn_simul ordn_simul',
       'simul month_simul', 'in_dur first_simul',
       'oneday_simul day_simul', 'month_simul simul',
       'cal_simul simul', 'cal_simul cal_simul', 
       'simul cal_simul', 'in_dur year_simul') 
    def multi_simul(self, p):
        getattr(p, 'in_dur', {})['function'] = 'telic_ext'
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

    @_('antdur_simul simul', 'antdur_simul year_simul')
    def multi_simul(self, p):
        p[0]['function'] = 'simultaneous'
        data = {
            'function': 'multi_simul',
            'parts': [p[0], p[1]],
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

    # (cl)
    @_('anterior_dur atelic_ext')
    def antdur_dur(self, p):
        data = {
            'function': 'antdur_dur',
            'parts': list(p),
        }

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

    # very complicated phrase!
    @_('day_simul simul begin_to_end',
       'hab_simul anterior_dist')
    def multi(self, p):
        data = {
            'function': '?',
            'parts': list(p),
        }
        return data

    # located durations
    @_('atelic_ext simul', 'atelic_ext year_simul')
    def dur_simul(self, p):
        data = {
            'function': 'dur_simul',
            'parts': [p[0], p[1]]
        }
        return data

    # (cl)
    @_('year_simul atelic_ext', 'day_simul atelic_ext',
       'in_dur atelic_ext', 'simul atelic_ext',
       'cal_simul multi_begintoend', 'cal_simul atelic_ext')
    def simul_dur(self, p):
        getattr(p,'in_dur',{})['function'] = 'simultaneous'
        data = {
            'function': 'simul_dur',
            'parts': list(p),
        }
        return data

    # (cl)
    @_('posterior atelic_ext')
    def postdur_dist(self, p):
        data = {
            'function': 'postdur_dist',
            'parts': list(p),
        }
        return data

    # simul phrases followed by potential references
    @_('simul antdur_simul', 'in_dur antdur_simul')
    def simul_ref(self, p):
        data = {
            'function': 'simul_ref',
            'parts': list(p),
        }
        return data

    # -- atelic extent
    @_('Ø duration', 'Ø num_year', 'Ø num_day', 
       'Ø one_year')
    def atelic_ext(self, p):
        getattr(p, 'one_year', {})['time'] = 'duration'
        p[1].update({
            'function': 'atelic_ext',
        })
        return p[1]

    # NB: this pattern is very significant
    # as it corroborates Haspelmath's hypothesis
    # that the zero-marked atelic extent derives
    # from the accusative / object role
    @_('>T duration', '>T num_day')
    def atelic_ext(self, p):
        getattr(p, 'num_day', {})['time'] = 'duration'
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
       '<L time', 'B person', 'Ø ordinal')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
        })
        return p[1] 

    @_('<L person')
    def simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': [
                {'time': 'duration', 'reference': 'person'},
            ],
        }
        return data

    @_('B simul', 'K simul', 'L simul',
       'K first_simul')
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

    @_('BEGINNING duration', 'BEGINNING time',
       'BEGINNING person',)
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
 
    @_('END duration', 'END time', 
       'END num_year', 'END num_day',
       'END person')
    def simul(self, p):
        p[1].update({
            'function': 'simultaneous',
            'time': 'singular',
            'time_loc': 'end',
            'count': 'measurement',
        })
        return p[1]

    # -- CALENDRICAL TIMES -- 
    @_('B month')
    def month_simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': [p[1]]
        }
        return data

    @_('B year', 'B num_year')
    def year_simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': [p[1]],
        }
        return data

    @_('year_simul num_year')
    def year_simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': list(p),
        }
        return data

    @_('B day', 'B num_day', 'B CARD', 
       'B THE CARD', 'Ø CARD')
    def day_simul(self, p):
        data = {
            'function': 'simultaneous',
            'time': 'single',
            'count': 'calendrical',
        }
        return data

    @_('day_simul year_simul')
    def day_simul(self, p):
        p[1]['time'] = 'duration'
        data = {
            'function': 'simultaneous',
            'parts': [
                {'time': 'single', 'count': 'calendrical'},
                p[1]
            ]
        }
        return data

    @_('B one_day')
    def oneday_simul(self, p):
        p[1].update({
            'function': 'simultaneous',
        })
        return p[1]

    @_('B ordinal')
    def ordn_simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': [p[1]]
        }
        return data

    @_('B first')
    def first_simul(self, p):
        data = {
            'function': 'simultaneous',
            'parts': [p[1]],
        }
        return data

    # -- CALENDRICAL COMPOSITIONS --
    @_('day_simul month_simul', 
       'day_simul ordn_simul year_simul',
       'month_simul day_simul', 
       'month_simul ordn_simul year_simul',
       'month_simul year_simul',
       'month_simul year_simul day_simul',
       'month_simul year',
       'ordn_simul day_simul',
       'first_simul day_simul',
       'year_simul month_simul day_simul',
       'year_simul ordn_simul day_simul',
       'year_simul month_simul', 
        # NB: last month_simul below belongs to day_simul
        # this is a work-around, but might be a better model
       'year_simul month_simul day_simul month_simul',
       'year_simul month_simul oneday_simul',
       'year_simul first_simul day_simul',
      )
    def cal_simul(self, p):
        getattr(p, 'ordn_simul', {})['number'] =  'calendrical'
        getattr(p, 'first_simul', {})['number'] =  'calendrical'
        getattr(p, 'year', {})['function'] =  'simultaneous'
        data = {
            'function': 'simultaneous',
            'parts': list(p),
            'reference': p[-1], # TODO: fix this properly
        } 
        return data

    # -- either atelic ext or simul --
    @_('Ø dur_sing', 'Ø one_day', 'Ø year',
       'Ø month')
    def atelic_simul(self, p):
        p[1].update({
            'function': 'atelic_ext, simultaneous',
        })
        return p[1]
    
    # -- telic extent / distances --
    @_('B duration', 'B dur_sing', 'K duration', 'B one_year')
    def in_dur(self, p):
        getattr(p, 'one_year', {})['time'] = 'duration'
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
    @_('atelic_ext anterior')
    def anterior_dist(self, p):
        data = {
            'function': 'anterior_dist',
            'parts': list(p), 
        }
        return data

    # -- anterior durative --
    @_('<D duration', '<D time', '<D simul',
       'L posterior', '<D dur_sing', '<D year',
       '<D num_year', '<D num_day', '<D one_day',
       '<D month', '<D person')
    def anterior_dur(self, p):
        getattr(p, 'num_year', {})['time'] = 'duration'
        getattr(p, 'num_day', {})['time'] = 'duration'
        data = {
            'function': 'anterior_dur',
            'parts': [p[1]],
        }
        return data

    @_('antdur_simul anterior_dur')
    def anterior_dur(self, p):
        p[0]['function'] = 'anterior_dur'
        p[0]['parts'] = [p[0]['time'], p[1]]
        return p[0]

    # -- anterior dur / simultaneous
    @_('L time', 'L duration', '>L time', 
       'L first', 'L num_day', 'L one_day')
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

    @_('>XR/ PERSON')
    def posterior(self, p):
        data = {
            'function': 'posterior',
            'reference': 'person',
        }
        return data

    # stand-alone >XR, will only match if 
    # followed by nothing else; causes shift/reduce
    # conflict, but this is tolerable for this proj.
    @_('>XR/')
    def posterior(self, p):
        return {
            'function': 'posterior',
            'advb': True,
        }

    @_('atelic_ext posterior')
    def posterior_dist(self, p):
        data = {
            'function': 'posterior_dist',
            'parts': list(p), 
        }
        return data

    # -- posteriors (durative?) --
    @_('MN time', 'MN dur_sing', 'MN one_day',
       'MN first', 'MN year', 'MN num_year', 'MN month')
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

    @_('MN simul', 'MN first_simul')
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

    @_('NUM duration', 'NUM time',
       'NUM_ONE duration')
    def duration(self, p):
        p[1].update({
            'time': 'durative',
            'count': 'mensural',
        })
        return p[1]

    @_('ALL duration', 'ALL time',
       'MANY duration', 'MANY time',
       'ALL year')
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

    @_('num_year num_year', 'duration num_year', 
       'num_day num_day', 'num_day duration',
       'duration num_day', 'duration year',
       'duration month')
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

    @_('GENDUR time', 'GENDUR month', 'GENDUR duration')
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
    
    # 2 Chr 21:15, ימים על־ימים
    # a kind of 'intensive' adjective
    @_('duration <L duration')
    def duration(self, p):
        data = {
            'time': 'durative',
            'parts': [p[0], p[2]],
            'quality': 'intensive',
        }
        return data

    @_('year >XR/ year')
    def duration(self, p):
        data = {
            'time': 'durative',
            'parts': [p[0], p[2]],
            'quality': 'intensive',
        }
        return data

    @_('NUM_ONE time')
    def dur_sing(self, p):
        p[1].update({
            'time': 'durative, singular', # disambig needed  
            'count': 'mensural',
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

    # define specialized time
    @_('BEGINNING', 'END')
    def time(self, p):
        return {
            'time': 'single',
            #'time_loc': 'beginning'
        }

    # calendrical time references
    @_('DAY')
    def day(self, p):
        return {
            'time': 'single',
        }

    @_('THE day')
    def day(self, p):
        p[1].update({
            'reference': 'anaphora',
            'ref_type': 'exophoric',
        })
        return p[1]

    @_('NUM_ONE day')
    def one_day(self, p):
        return {
            'time': 'single',
        }

    @_('NUM day')
    def num_day(self, p):
        p[1].update({
            'number': 'calendrical',
        })
        return p[1]

    @_('num_day month_ref', 'num_day person_ref',
       'one_day month_ref', 'day month_ref')
    def day(self, p):
        data = {
            'parts': [p[0], p[1]],
            'reference': p[1]['reference'],
        }
        return data

    @_('CARD month_ref', 'CARD ordn_ref')
    def day(self, p):
        getattr(p, 'ordn_ref', {})['reference'] = 'month'
        data = {
            'parts': [
                {'time': 'single', 'count': 'calendrical'},
                p[1]
            ],
            'reference': p[1]['reference'],
        }
        return data

    @_('MONTH')
    def month(self, p):
        return {
            'time': 'single'
        }

    @_('THE month')
    def month(self, p):
        p[1].update({
            'reference': 'anaphora',
            'ref_type': 'exophoric',
        })
        return p[1]

    @_('THIS month')
    def month(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'near'
        })
        return p[1]

    @_('THAT month')
    def month(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'far'
        })
        return p[1]

    @_('NUM month')
    def month(self, p):
        p[1].update({
            'number': 'calendrical',
        })
        return p[1]

    @_('SFX month')
    def month(self, p):
        p[1].update({
            'reference': 'personal',
        })
        return p[1]

    @_('month antdur_simul')
    def month(self, p):
        data = {
            'reference': 'time',
            'parts': [
                {'time': p[1]['time']}
            ]
        }
        return data

    @_('month L year')
    def month(self, p):
        data = {
            'reference': 'year',
            'parts': [
                p[0], p[2]
            ]
        }
        return data

    @_('GENCARD day')
    def day(self, p):
        p[1].update({
            'time': 'single',
            'count': 'calendrical',
        })
        return p[1]

    @_('YEAR')
    def year(self, p):
        return {
            'time': 'single'
        }

    @_('THE year')
    def year(self, p):
        p[1].update({
            'reference': 'anaphora',
            'ref_type': 'exophoric',
        })
        return p[1]

    @_('THIS year')
    def year(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'near'
        })
        return p[1]

    @_('THAT year')
    def year(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'spatial',
            'ref_dist': 'far'
        })
        return p[1]

    @_('GENCARD year')
    def year(self, p):
        p[1].update({
            'time': 'single',
            'count': 'calendrical',
        })
        return p[1]

    @_('NUM year')
    def num_year(self, p):
        p[1].update({
            'number': 'calendrical',
        })
        return p[1]

    @_('NUM_ONE year')
    def one_year(self, p):
        return {
            'time': 'single',
        }

    @_('year person_ref', 'num_year person_ref')
    def year(self, p):
        data = {
            'parts': [p[1]],
            'reference': 'person',    
        }
        return data
    
    @_('year antdur_simul')
    def year(self, p):
        data = {
            'reference': 'time',
            'parts': [
                {'time': p[1]['time']}
            ]
        }
        return data

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

    @_('ORDN')
    def ordinal(self, p):
        return {
            'time': 'singular',
            'reference': 'absolute', # ??
        }

    @_('THE ordinal')
    def ordinal(self, p):
        return p[1]

    # see Gen 28:19 for why this pattern is needed:
    # לראשנה
    # to treat RISHON as a normal ordinal would mean
    # that this phrase would get interpreted as a 
    # calendrical reference; whereas it never occurs
    # as such in calendrical simultaneous phrases
    @_('FIRST')
    def first(self, p):
        return {
            'time': 'singular',
            'reference': 'absolute',
        }
    @_('THE first')
    def first(self, p):
        return p[1]

    @_('PERSON')
    def person(self, p):
        return {
            'reference': 'person',
        }
    @_('THE person')
    def person(self, p):
        p[1].update({
            'reference': 'anaphora',
            'ref_type': 'exophoric',
        })
        return p[1]
    @_('SFX person')
    def person(self, p):
        p[1].update({
            'reference': 'deictic',
            'ref_type': 'personal',
        })
        return p[1]
    
    # define reference constructions
    @_('L person')
    def person_ref(self, p):
        return {
            'reference': 'person', # TODO: convert to list to keep person refs
        }

    @_('L month',
       'IT B')
    def month_ref(self, p):
        return {
            'time': 'single',
            'reference': 'month'
        }
    
    @_('L ordinal')
    def ordn_ref(self, p):
        return {
            'time': 'single',
            'reference': 'absolute',
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

    F, L = API.F, API.L
    
    # gather eligible clauses
    # an eligible clause is one for which all
    # its phrases are parsed
    clauses = {}
    good_phrases = set(phrases)
    for clause in F.otype.s('clause'):
        time_phrases = [
            ph for ph in L.d(clause, 'phrase')
                if F.function.v(ph) == 'Time'
        ]
        if time_phrases and good_phrases.issuperset(set(time_phrases)):
            clauses[clause] = time_phrases

    print(len(clauses), 'clauses ready for parsing...')

    for clause, time_phrases in clauses.items(): 

        # collect tokens for all phrase parsings
        tokens = []
        for ph in time_phrases:
            parsing = phrases[ph]
            tokens.extend(
                tokenizer.tokenize(parsing)
            )

        # attempt parsing
        parsing = parser.parse(t for t in tokens)
        if parsing is not None and error_tracker['e'] is None:
            parsed[clause] = parsing
        else:
            toks = str([t.type for t in tokens])
            errors[clause] = (error_tracker['e'], toks)

            # reset error tracker
            error_tracker['e'] = None 

    # export
    with open(paths['parsed'], 'w') as outfile:
        json.dump(parsed, outfile, indent=2, ensure_ascii=False)

    with open(paths['notparsed'], 'w') as outfile:
        json.dump(errors, outfile, indent=2, ensure_ascii=False)
