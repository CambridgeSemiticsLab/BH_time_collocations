import json
from sly import Parser as SlyParser
from sly.lex import Token as SlyToken
from data.pipeline.legacy_scripts.tools import nav_tree as nt
from data.pipeline.legacy_scripts.tools.load_parse import ParseLoader


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

        # set some linguistic standards

        self.seqmap = {
            '>XRWN/': 'POST',
            '>XR=/': 'POST',
            'TJKWN/': 'MIDDLE',
        }
        self.tensemap = {
            'MXR/': 'FUT',
        }
        self.demon_map = {
            "Z>T": "THIS",
            "HJ>": "THAT",
            "HMH": "THAT",
            ">LH": "THIS",
            "HM": "THAT",
            "HW>": "THAT",
            "ZH": "THIS",
        }
    
    def sly_token(self, tag, parse, **values):
        """Get SLY token object with customized source_data."""
        token = SlyToken()
        self.tokens.add(tag) # track unique tokens
        token.value = dict(values)
        token.value['parse'] = parse
        slots = sorted(self.get_slots(parse))
        token.value['phatom'] = self.L.u(slots[0], 'phrase_atom')[0]
        token.value['slots'] = slots
        token.type = tag
        token.index = self.i
        token.lineno = 1
        self.i += 1

#        if token.value['phatom'] in {906675, 906676, 906677, 906678, 906679, 906680}:
#            print(f'phatom={token.value["phatom"]}; val={tag}, index={token.index}')
#            print(f'\t{slots}')

        return token

    def get_head(self, phrase):
        """Error-free head grabbing, allowing for singular phrases."""
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

    def tokenize(self, parsedphrase, start=True, nulltoken=False, prefixes=[]):
        """Follow path to right-most item (head) and yield tokens.
        
        tokenize must adjudicate what gets yielded as a token
        and what does not. At times, we need to direct tokenize down
        a path that is not on the head-path (e.g. parallel relas).
        Sometimes whether a rela becomes a token is conditional 
        on a number of factors tuned to produce good tokens.
        """
        
        # ignore these relations
        ignore = {None, 'ADVB',}
        def tag_rela(rela, name):
            return rela == name and rela not in ignore
   
        # process all internal relations
        head_phrases = self.get_head_phrases(parsedphrase)
#        if self.debug and start:
#            print('START')
#            pprint(head_phrases)
#            print()
        relas = set(self.get_rela(ph)[-1] for ph in head_phrases)
        paras = []

        # reset token indexer
        if start:
            self.i = 0

        # yield any before-phrase tokens
        # e.g. TIMEAPPO
        for pre in prefixes:
            yield pre
        
        first_sp = True
        for phrase in head_phrases:

            src, tgt, rela = self.get_rela(phrase) 

            # yield null token if no starting preposition found
            if first_sp and (start or nulltoken):
                token = self.null_token(phrase)
                if token:
                    yield token

            first_sp = False

            # tokenize prepositions
            if tag_rela(rela, 'PP'):
                yield self.prep_token(phrase)
                
            # handle parallels
            elif tag_rela(rela, 'PARA'):
                for subphrase in nt.unfold_paras(src):
                    paras.extend(
                        self.tokenize(subphrase, start=False) # recursively tokenize them
                    )                    
            
            # handle phrases which should be parsed in total
            elif tag_rela(rela, 'SPEC'):
                if type(src) == int:
                    src = [None, src, None]
                paras.extend(
                    self.tokenize(
                        src, 
                        start=False,
                    )
                )
    
            elif tag_rela(rela, 'TIMEAPPO'):
                token = self.timeappo_token(phrase)
                paras.extend(
                    self.tokenize(
                        src,
                        start=False,
                        nulltoken=True,
                        prefixes=[token]
                    )
                )

            # tokenize appositional relations if relevant
            elif tag_rela(rela, 'APPO'):
                token = self.appo_token(phrase, paras)
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

            elif tag_rela(rela, 'ADJV'):
                token = self.adjv_token(phrase)
                if token:
                    yield token

            # genitival relations
            elif tag_rela(rela, 'GP'):
                token = self.gen_token(phrase, paras)
                if token:
                    yield token

            elif tag_rela(rela, 'CARDC'):
                yield self.sly_token(rela, phrase)

           # TODO: test if this can now be eliminated safely
           # drip-bucket tokenizer
            elif rela not in ignore:
                yield self.sly_token(rela, src)

            # give TIME token once reaching the end
            if head_phrases[-1] == phrase:
                timetoks = list(self.time_token(phrase, relas))
                seenlex = set(tok.value.get('lex') for tok in timetoks)
                for tok in timetoks:
                    yield tok

                # yield all of the parallel tokens in sorted order
                sort_paras = lambda p: [p.value['phatom'], p.index]
                sorted_paras = sorted(paras, key=sort_paras)
#                if self.debug:
#                    print([p.type for p in paras])

                for sp_token in sorted_paras:

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
                yield self.sly_token('IT', time) # בו
            else:
                yield self.sly_token('SFX', time)
        if F.nu.v(time) == 'du':
            yield self.sly_token('NUM', time)
        if F.uvf.v(time) == 'H':
            yield self.sly_token('LOCALE', time)

        # get source_data for processing time words
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
        yield self.sly_token(token, time, lex=lex)
    
    def prep_token(self, phrase):
        """Parse prepositions into singular tokens."""
        token = self.F.lex.v(phrase[0])
        token = self.lexmap.get(token, token)
        return self.sly_token(token, phrase[0])
    
    def null_token(self, phrase):
        """Parse whether phrase begins without preposition."""
      
        # process single-word phrases
        if len(phrase) == 1:
            if self.slot2pos[phrase[0]] != 'PREP':
                return self.sly_token('Ø', phrase[0])
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
            return self.sly_token('Ø', slots[0])

    def def_token(self, phrase):
        """Definite tokens."""
        return self.sly_token('THE', phrase[0])

    def appo_token(self, phrase, paras):
        """Parse appositional tokens."""
        src, tgt, rela = self.get_rela(phrase)
        
        # get the head item of the appositional phrase
        appo_head = self.get_head(src)
        appo_lex = self.F.lex.v(appo_head)
        head = self.get_head(tgt)
        head_lex = self.F.lex.v(head)
       
        # process appositional demonstratives
        if self.F.pdp.v(appo_head) == 'prde':
            token = self.demon_map[appo_lex]
            return self.sly_token(token, appo_head)

        elif self.slot2pos[appo_head] == 'ORDN':
            return self.sly_token('ORDNM', appo_head) 

        elif appo_lex == 'RB/':
            return self.sly_token('MANY', appo_head)

        elif appo_lex in self.seqmap:
            # (name of value to modified, value itself)
            mod = ('seqs', self.seqmap[appo_lex])
            return self.sly_token('MOD', appo_head, mod=mod)

        elif appo_lex in self.tensemap:
            mod = ('tenses', self.tensemap[appo_lex])
            return self.sly_token('MOD', appo_head, mod=mod)

        elif appo_lex == 'CLCWM/':
            return self.sly_token('DISTAGO', appo_head)

        elif (
            appo_lex == head_lex
            and type(src) == int
            and type(tgt) == int
        ):
            paras.extend(
                self.time_token(
                    [None, src, None],
                    {}
                )
            )

    def demon_token(self, phrase):
        """Parse demonstratives."""
        src, tgt, rela = self.get_rela(phrase)
        demon = self.get_head(src)
        demonlex = self.F.lex.v(demon)
        token = self.demon_map[demonlex]
        return self.sly_token(token, demon)

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
        return self.sly_token(token, phrase[0])

    def adjv_token(self, phrase):
        """Parse adjectival phrase tokens."""
        src, tgt, rel = self.get_rela(phrase)
        adjv_head = self.get_head(src)
        adjv_lex = self.F.lex.v(adjv_head)
        tgt_head = self.get_head(tgt)
    
        if (
            type(tgt) == list 
            and tgt[-1] == 'NUM'
            and adjv_lex == '<WD/'
        ):
            return self.sly_token('YETDUR', src)    
        elif (
            type(src) == list
            and type(src[0]) == int
            and self.slot2pos[adjv_head] in {'CARD', 'CARD1'}
            and self.F.lex.v(src[0]) == '<WD/'
        ):
            return self.sly_token('YETDUR', src)
        elif (
            adjv_lex in self.seqmap 
        ):
            mod = ('seqs', self.seqmap[adjv_lex])
            return self.sly_token('MOD', src, mod=mod)
        elif (
            adjv_lex in self.tensemap
        ):
            mod = ('tenses', self.tensemap[adjv_lex])
            return self.sly_token('MOD', src, mod=mod)
        elif(
            adjv_lex in self.demon_map
        ):
            tok = self.demon_map[adjv_lex]
            return self.sly_token(tok, src)

    def gen_token(self, phrase, paras):
        """Parse genitive phrase tokens."""
        src, tgt, rel = self.get_rela(phrase)
        gen_head = self.get_head(src)
        ph_head = self.get_head(tgt)
        F = self.F

        # identify genitive durations
        # e.g. חדשׁ ימים "month of days"
        dur_lexs = {'JWM/', 'DWR/', 'NYX/'}

        if (
            self.F.nu.v(ph_head) == 'sg'
            and self.F.lex.v(gen_head) in dur_lexs
            and self.F.nu.v(gen_head) == 'pl'
        ):
            return self.sly_token('GENDUR', src)

        elif (
            self.slot2pos[gen_head] in {'CARD1', 'CARD'}
        ):
            return self.sly_token('GENCARD', src)

        elif(
            type(src) == list
            and src[-1] == 'NUM'
            and F.lex.v(gen_head) == F.lex.v(ph_head) 
        ):
            return self.sly_token('NUM', src)
        elif(
            self.slot2pos[gen_head] == 'ORDN'
        ):
            return self.sly_token('ORDNM', src)
        elif(
            F.nametype.v(gen_head) == 'pers'
            and F.nametype.v(ph_head) != 'pers'
            and F.lex.v(ph_head) not in {'MLKWT/'}
        ):
            return self.sly_token('PERSONM', src)
        elif (
            F.lex.v(gen_head) == 'XDC=/'
            and F.lex.v(ph_head) == 'JWM/'
        ):
            paras.append(self.sly_token('MONTHREF', src))
            paras.extend(self.tokenize(src, start=False))
        else:
            return self.sly_token('GENREF', src) 

    def quant_token(self, phrase):
        """Return quantifier tokens."""
        src, tgt, rel = self.get_rela(phrase)
        quant = self.get_head(src)
        quantlex = self.F.lex.v(quant)
        token = self.lexmap.get(quantlex)
        if token:
            return self.sly_token(token, src)

    def timeappo_token(self, phrase):
        """Return time appositional token."""
        src, tgt, rel = self.get_rela(phrase)
        return self.sly_token('TIMEAPPO', src)

# TimeParser helper functions for organizing the source_data

def expand(time, key, val):
    """Append to a list-based value in the time dicts."""
    try:
        time[key].insert(0, val)
    except KeyError:
        time[key] = [val]
    except TypeError:
        raise Exception(time, key, val)

def expand_key(key):
    """Return func for expanding a given key."""
    def key_adder(val_str, val, time):
        """Add source_data to a time based on key/val."""

        # determine if the item is a parsed time or a token
        if 'times' in val:
            obj = val
        else: 
            obj = val['parse']

        # append the tag to the appropriate key, its val is a list
        tag = [val_str, obj]
        expand(time, key, tag) 
        time['slots'].extend(val['slots'])
        time['slots'] = sorted(set(time['slots']))
        return time

    return key_adder

add_ref = expand_key('refs')
add_quant = expand_key('quants')
add_prep = expand_key('preps')
add_seq = expand_key('seqs')
add_tense = expand_key('tenses')

def add_function(funct, time):
    """Add a new function."""
    expand(time, 'functions', funct)
    return time

def add_qual(qual, time):
    """Add new temporal qualities."""
    expand(time, 'quals', qual)
    return time

def init_time(p, **kwargs):
    """Initialize time source_data."""
    time_data = {
        'times': [p['parse']],
        'slots': p['slots'],
    }    
    time_data.update(**kwargs)
    return time_data

def init_times(times, **kwargs):
    """Add attributes to a set of times."""
    for t in times:
        t.update(**kwargs)
    return times

def init_ctime(*times, **kwargs):
    """Initialize new time composed of times."""
    time_data = {
        'times': [],
        'slots': [],
        'functions': [],
    }
    for time in times:
        time_data['times'].append(time)
        time_data['slots'].extend(time['slots'])
#        time_data['functions'].extend(time.get('functions', []))
    time_data.update(**kwargs)
    return time_data

def mergetime(small, big):
    """Merge one time into another."""
    big['slots'].extend(small['slots'])
    big['times'].append(small)
    return big

def conv_num(numtype, time, as_dur=False):
    """Convert num to a calendar num or quant."""
    if numtype == 'CALNUM':
        add_ref(numtype, time['NUM'], time)
    elif numtype == 'NUMQ':
        add_quant(numtype, time['NUM'], time)
        if as_dur:
            add_qual('durative', time)

    del time['NUM']
    return time

def conv_nums(numtype, times, **kwargs):
    """Convert multiple nums at once."""
    for time in times:
        if type(time) == dict and 'NUM' in time:
            conv_num(numtype, time, **kwargs)

def check_nums(time):
    """Checks for conditions to convert numbers."""
    refs = set(r[0] for r in time.get('refs', []))
    if {'PERS', 'YEAR', 'MONTH'} & refs:
        conv_num('CALNUM', time)

def add_mod(modtoken, time):
    """Add modifier to a time."""
    key, val_str = modtoken['mod']
    del modtoken['mod']
    if key == 'tenses':
        add_tense(val_str, modtoken, time)
    elif key == 'seqs':
        add_seq(val_str, modtoken, time)
    else:
        raise Exception(f'modifier key {key} not yet added!')
    return time

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
        '<M',
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
        'OFTEN',
        'THUS',
        'THEN',
        'CARD',
        'ONWARD',
        'MONTH',
        'YEAR',
        'DAY',
        'PERSON',
        'PERSONM',
        'IT',
        'LOCALE',
        'ORDN',
        'ORDNM',
        'FIRST',
        'MOD',
        'GENREF',
        'TIMEAPPO',
        'YETDUR',
        'DISTAGO',
        'MONTHREF',
        'CARDC'
    }

    def error(self, token):
        """Keep track of errors."""
        try:
            self.error_tracker['e'] = token.type
        except:
            self.error_tracker['e'] = 'reached end'

    debugfile = '../results/data_metrics/timeparse.debug'
       
    # -- FINAL MATCHES --
    @_('timephrase', 'timeappo')
    def final(self, p):
        return p[0]

    @_('atelic_ext', 'simul', 'in_dur',
       'anterior', 'anterior_dur', 'posts',
       'posterior', 'atelic_simul', 'antdur_simul',
       'hab_simul', 'begin_to_end', 'multi_simul',
       'antpost_dur', 'posterior_dur', 'multi_postantdur',
       'anterior_posterior', 'multi_antdursim', 'dur_to_end',
       'simul_posts', 'multi_posts', 'simul_to_end', 
       'dur_simul', 'month_ref', 
       'cal_simul', 'habitual', 'multi_habitual',
       'multi',  'simul_ref', 'dist_prosp',
       'multi_begintoend', 'posterior_dist', 'anterior_dist',
       'simul_dur', 'postdur_dist', 'antdur_dur', 'dist_fut',
       'distfut_dur', 'dist_past', 'multi_simuls', 
       'simuls_or_ext', 'multi_antdur')
    def timephrase(self, p): 
        return p[0] 

    @_('month_simul', 'year_simul','day_simul',)
    def timephrase(self, p):
        conv_nums('CALNUM', p)
        return p[0]

    @_('oneday_simul')
    def timephrase(self, p):
        conv_nums('NUMQ', p)
        return p[0]

    @_('ordn_simul')
    def timephrase(self, p):
        add_ref('ORDN', p[0]['ordn'], p[0])
        del p[0]['ordn']
        return p[0]

    @_('first_simul')
    def timephrase(self, p):
        del p[0]['ordn'] # treat it as stand-alone word
        return p[0]

    @_('timephrase TIMEAPPO timephrase')
    def timeappo(self, p):
        p[0]['slots'].extend(p[2]['slots'])
        p[0]['timeappos'] = [p[2]]
        return p[0]

    @_('timephrase TIMEAPPO timeappo')
    def timeappo(self, p):
        p[0]['slots'].extend(p[2]['slots'])
        p[0]['timeappos'] = [p[2]]
        return p[0]
         
    # -- composed constructions
    @_('simul simul')
    def hab_simul(self, p):
        return init_ctime(
            *p,    
            functions=['habitual, multi_simul'],
            quals=['distributive'],
        )
        
    @_('hab_simul in_dur') 
    def hab_simul(self, p):
        return mergetime(p[1], p[0])

    @_('simul hab_simul')
    def hab_simul(self, p):
        return mergetime(p[0], p[1])

    @_('hab_simul hab_simul')
    def hab_simul(self, p):
        return mergetime(p[1], p[0])

    @_('Ø card antdur_simul')
    def hab_simul(self, p):
        p[2]['functions'][0] = 'simultaneous'
        add_function('simultaneous', p[1])
        return init_ctime(
            p[1], p[2],
            functions=['habitual, multi_simul'],
            quals=['distributive'],
        )    
    
    @_('Ø year year_simul', 'Ø day day_simul', 'Ø month month_simul')
    def habitual(self, p):
        add_function('simultaneous', p[1])
        return init_ctime(
            p[1], p[2],
            functions=['habitual'],
            quals=['distributive'],
        )

    @_('MN OFTEN year year_simul', 'MN OFTEN month month_simul',
       'MN OFTEN time simul')
    def habitual(self, p):
        time = init_ctime(
            p[2], p[3],
            functions=['habitual'],
            quals=['distributive'],
        )
        add_quant('OFTEN', p[1], time)
        add_prep('MN', p[0], time)
        return time

    @_('habitual habitual')
    def multi_habitual(self, p):
        return init_ctime(
            p[0], p[1],
            functions=['multi_habitual']
        )

    @_('posterior_dur anterior_dur', 'antpost_dur anterior_dur')
    def begin_to_end(self, p):
        p[0]['functions'][0] = 'posterior_dur'
        return init_ctime(
            *p,
            functions=['begin_to_end'],
            quals=['durative'],
        )

    @_('posts anterior_dur', 'posts antdur_simul')
    def begin_to_end(self, p):
        p[0]['functions'][0] = 'posterior_dur'        
        add_qual('durative', p[0])
        p[1]['functions'][0] = 'anterior_dur'
        return init_ctime(
            *p,
            functions=['begin_to_end'],
            quals=['durative'],
        )

    @_('posts ONWARD')
    def begin_to_end(self, p):
        p[0]['functions'][0] = 'posterior_dur'
        add_qual('durative', p[0])
        onward = init_time(
            p[1], 
            quals=['durative'],
            functions=['anterior_dur']
        )
        return init_ctime(
            p[0], onward,
            functions=['begin_to_end'],
            quals=['durative'],
        )

    @_('posts LOCALE ONWARD')
    def begin_to_end(self, p):
        p[0]['functions'][0] = 'posterior_dur'
        add_qual('durative', p[0])
        onward = init_time(
            p[2], 
            quals=['durative'], 
            functions=['anterior_dur'],
            locale=True
        )
        return init_ctime(
            p[0], onward,
            functions=['begin_to_end'],
            quals=['durative'],
        )

    @_('posts LOCALE duration')
    def begin_to_end(self, p):
        p[0]['functions'][0] = 'posterior_dur'
        add_qual('durative', p[0])
        add_function('anterior_dur', p[2])
        p[2]['locale'] = True
        return init_ctime(
            p[0], p[2],
            functions=['begin_to_end'],
            quals=['durative'],
        )

    @_('MN month month_ref')
    def begin_to_end(self, p):
        add_prep('MN', p[0], p[1])
        add_function('posterior_dur', p[1])
        add_qual('durative', p[1])
        p[2]['functions'][0] = 'anterior_dur'
        return init_ctime(
            p[1], p[2],
            functions=['begin_to_end'],
            quals=['durative'],
        )

    @_('begin_to_end begin_to_end')
    def multi_begintoend(self, p):
        return init_ctime(
            *p,
            functions=['multi_begin_to_end'],
            quals=['durative'],
        )

    @_('atelic_ext anterior_dur')
    def dur_to_end(self, p):
        return init_ctime(
            *p,
            functions=['dur_to_end'],
            quals=['durative'],
        )

    @_('atelic_ext antdur_simul')
    def dur_to_end(self, p):
        p[1]['functions'][0] = 'anterior_dur'
        return init_ctime(
            *p,
            functions=['dur_to_end'],
            quals=['durative'],
        )

    # 'on the eigth day and onward'
    @_('in_dur anterior_dur', 'simul anterior_dur')
    def simul_to_end(self, p):
        return init_ctime(
            *p,
            functions=['simul_to_end'],
            quals=['durative'],
        )

    @_('simul ONWARD')
    def simul_to_end(self, p):
        onward = init_time(p[1], quals=['durative'])
        return init_ctime(
            p[0], onward,
            functions=['simul_to_end'],
            quals=['durative'],
        )

    @_('in_dur simul', 'simul in_dur', 'in_dur in_dur',
       'day_simul day_simul', 'day_simul simul',
       'year_simul simul', 'ordn_simul ordn_simul',
       'simul month_simul', 'simul first_simul',
       'oneday_simul day_simul', 'month_simul simul',
       'cal_simul cal_simul',
       'in_dur year_simul',
       'first_simul simul',
       'simul year_simul', 'simul oneday_simul',
    ) 
    def multi_simul(self, p):
        return init_ctime(
            *p,
            functions=['multi_simul'],
        )

    @_('cal_simul simul', 'simul cal_simul',
       'in_dur cal_simul')
    def multi_simul(self, p):
        conv_nums('CALNUM', p)
        return init_ctime(
            *p,
            functions=['multi_simul'],
        )

    @_('multi_simul simul') 
    def multi_simul(self, p):
        return mergetime(p[1], p[0])

    @_('simul multi_simul')
    def multi_simul(self, p):
        return mergetime(p[0], p[1]) 

    @_('antdur_simul simul', 'antdur_simul year_simul')
    def multi_simul(self, p):
        p[0]['functions'][0] = 'simultaneous'
        return init_ctime(
            *p,
            functions=['multi_simul'],
        )

    @_('Ø time time num_dur')
    def multi_simul(self, p):
        times = init_times(
            list(p)[1:],
            functions=['simultaneous'],
        )
        return init_ctime(
            *times,
            functions=['multi_simul'],
        )

    @_('Ø THE multi_time')
    def multi_simul(self, p):
        add_ref('THE', p[1], p[2]['times'][0])
        add_function('multi_simuls', p[2])
        for time in p.multi_time['times']:
            add_function('simultaneous', time)
        return p[2]

    @_('posts anterior')
    def multi_postantdur(self, p):
        p[0]['functions'][0] = 'anterior_dur'
        return init_ctime(
            *p,
            functions=['anterior_dur + anterior'],
        )

    @_('anterior posterior')
    def anterior_posterior(self, p):
        return init_ctime(
            *p,
            functions=['anterior + posterior'],
        )

    @_('antdur_simul antdur_simul')
    def multi_antdursim(self, p):
        return init_ctime(
            *p,
            functions=['multi_simul, multi_antdur'],
        )
   
    @_('antdur_simul multi_antdursim')
    def multi_antdursim(self, p):
        return mergetime(p[0], p[1])

    # (cl)
    @_('anterior_dur atelic_ext')
    def antdur_dur(self, p):
        return init_ctime(
            *p,  
            functions=['anterior_dur + duration'],
        )

    @_('anterior_dur all_dur')
    def antdur_dur(self, p):
        add_function('atelic_ext', p[1])
        add_qual('durative', p[1])
        return init_ctime(
            *p,  
            functions=['anterior_dur + duration'],
        )

    # 'in that day and after tomorrow'
    @_('simul posts')
    def simul_posts(self, p):
        return init_ctime(
            *p,  
            functions=['simul_posts'],
        ) 

    @_('posts posts')
    def multi_posts(self, p):
        return init_ctime(
            *p,
            functions=['multi_posts'],
        )

    @_('multi_posts posts')
    def multi_posts(self, p):
        return mergetime(p[1], p[0])

    # complicated phrases to be disambiguated by hand
    @_('day_simul simul begin_to_end',
       'antdur_simul habitual',
       'posterior simul simul', 'posts simul',
       'posts begin_to_end', 'posterior_dur simul',
       'habitual begin_to_end',
       'hab_simul begin_to_end', 'in_dur multi',
       'posterior antdur_simul',
       'posterior year_simul',
       'antdur_simul begin_to_end',
       'atelic_ext anterior atelic_ext'
    )
    def multi(self, p):
        return init_ctime(
            *p,
            functions=['?'],
        )

    @_('Ø time distago') # NB: distago here cannot be modifier due to position
    def multi_simuls(self, p):
        times = list(p)[1:]
        for time in times:
            add_function('simultaneous', time)
        return init_ctime(
            *times,
            functions=['multi_simuls']
        )

    @_('year_simul month_simul day_simul simul')
    def multi_simul(self, p):
        cal_simul = init_ctime(
            *list(p)[:-1],
            functions=['cal_simul'],
        )
        return init_ctime(
            cal_simul, p[-1],
            functions=['multi_simul'],
        )

    # located durations
    @_('atelic_ext simul', 'atelic_ext year_simul')
    def dur_simul(self, p):
       return init_ctime(
            *p,
            functions=['atelic_ext + simultaneous'],
        )

    @_('atelic_ext in_dur')
    def dur_simul(self, p):
        p[1]['functions'][0] = 'simultaneous'
        return init_ctime(
            *p,
            functions=['atelic_ext + simultaneous'],
        )

    # (cl)
    @_('year_simul atelic_ext', 'day_simul atelic_ext',
       'simul atelic_ext', 'cal_simul multi_begintoend', 
       'multi_simul atelic_ext')
    def simul_dur(self, p):
        return init_ctime(
            *p,
            functions=['simultaneous + atelic_ext']
        )         

    @_('Ø THE time num_dur')
    def simul_dur(self, p):
        add_ref('THE', p.THE, p.time)
        add_function('simultaneous', p.time)
        add_function('atelic_ext', p.num_dur)
        return init_ctime(
            p.time, p.num_dur,
            functions=['simultaneous + atelic_ext'],
        )

    @_('L THE time num_dur')
    def simul_dur(self, p):
        add_ref('THE', p.THE, p.time)
        add_prep('L', p.L, p.time)
        add_function('simultaneous', p.time)
        add_function('atelic_ext', p.num_dur)
        return init_ctime(
            p.time, p.num_dur,
            functions=['simultaneous + atelic_ext'],
        )

    @_('in_dur atelic_ext')
    def simul_dur(self, p):
        p[0]['functions'][0] = 'simultaneous'
        return init_ctime(
            *p,
            functions=['simultaneous + atelic_ext']
        )         

    @_('dist_fut atelic_ext')
    def distfut_dur(self, p):
        return init_ctime(
            *p,
            functions=['distfut_ext'],
        )

    # (cl)
    @_('posterior atelic_ext')
    def postdur_dist(self, p):
        return init_ctime(
            *p,
            functions=['posterior + atelic_ext'],
        )

    @_('Ø duration upon_year')
    def dist_prosp(self, p):
        p[2]['functions'][0] = 'reference'
        add_function('dist_prospective', p[1])
        add_ref('YEAR', p[2], p[1])
        return p[1]

    # simul phrases followed by potential references
    @_('simul antdur_simul', 'dist_fut antdur_simul')
    def simul_ref(self, p):
        p[1]['functions'][0] = 'reference'
        add_ref('LREF', p[1], p[0])
        return p[0]

    # -- atelic extent
    @_('Ø duration', 'Ø num_year', 
       'Ø num_day', 'Ø one_year',
       'Ø num_dur', 'Ø all_dur')
    def atelic_ext(self, p):
        conv_nums('NUMQ', p, as_dur=True)
        add_function('atelic_ext', p[1])
        return p[1]

    @_('Ø YETDUR num_dur')
    def atelic_ext(self, p):
        conv_nums('NUMQ', p, as_dur=True)
        add_function('atelic_ext', p[2])
        return p[2]

    # ambiguous
    @_('Ø multi_time') 
    def simuls_or_ext(self, p):
        add_function('multi_simuls, atelic_ext', p[1])
        for time in p[1]['times']:
            add_function('simultaneous, atelic_ext', time)
        return p[1]

    # NB: this pattern is very significant
    # as it corroborates Haspelmath's hypothesis
    # that the zero-marked atelic extent derives
    # from the accusative / object role
    @_('>T num_dur', '>T num_day')
    def atelic_ext(self, p):
        conv_nums('NUMQ', p, as_dur=True)
        add_function('atelic_ext', p[1])
        add_prep('>T', p[0], p[1])
        return p[1]

    @_('<M time')
    def atelic_ext(self, p):
        add_function('atelic_ext', p[1])
        add_prep('<M', p[0], p[1])
        return p[1]

    @_('atelic_ext atelic_ext')
    def atelic_ext(self, p):
        time = init_ctime(
            *p,
            functions=['atelic_ext']
        )
        return time

    @_('Ø num_dur month duration', 'Ø time num_dur')
    def atelic_ext(self, p):
        times = init_times(
            list(p)[1:],
            functions=['atelic_ext'],
        )
        return init_ctime(
            *times,
            functions=['atelic_ext'], 
        )

    # -- future / past distances --
    @_('distago time')
    def dist_past(self, p):
        quants = p.distago['quants']
        tenses = p.time['tenses']
        return init_ctime(
            *p,
            functions=['dist_past'],
            quants=quants,
            tenses=tenses,
        ) 

    @_('Ø dist_past')
    def dist_past(self, p):
        return p[1]

    @_('B YETDUR num_dur', 'B YETDUR num_year')
    def dist_fut(self, p):
        conv_nums('NUMQ', p)
        add_tense('FUT', p.YETDUR, p[-1])
        add_prep('B', p.B, p[-1])
        add_function('dist_fut', p[-1])
        return p[-1]

    @_('L YETDUR duration')
    def dist_fut(self, p):
        add_tense('FUT', p.YETDUR, p[-1])
        add_prep('L', p.L, p[-1])
        add_function('dist_fut', p[-1])
        return p[-1]

    # -- simultaneous --
    @_('B time', 'B person', 'B duration', 'B all_dur')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1] 

    @_('B multi_time')
    def simul(self, p):
        add_qual('durative', p[1])
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1] 

    @_('BJN/ num_dur')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('BJN/', p[0], p[1])
        return p[1] 

    @_('K time', 'K duration', 'K num_dur', 'K dist_past')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('K', p[0], p[1])
        return p[1] 

    @_('<L time', '<L person')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('<L', p[0], p[1])
        return p[1] 

    @_('<L year')
    def upon_year(self, p):
        add_function('simultaneous', p[1])
        add_prep('<L', p[0], p[1])
        return p[1]

    @_('Ø time')
    def simul(self, p):
        add_function('simultaneous', p[1])
        return p[1] 

    @_('Ø ordinal')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_ref('ORDN', p[1]['ordn'], p[1])
        del p[1]['ordn']
        return p[1]

    @_('B simul')
    def simul(self, p):
        add_prep('B', p[0], p[1])
        add_function('simultaneous', p[1])
        return p[1]

    @_('K simul', 
       'K first_simul')
    def simul(self, p):
        add_prep('K', p[0], p[1])
        add_function('simultaneous', p[1])
        return p[1]

    # Very interesting pattern here, somewhat
    # similar to antpost_dur:
    # למן־היום (II Sam 19:25)
    # except here I think this is a simultaneous
    # position after a duration: 
    # כמשלש חדשים (Gen 38:24)
    @_('K posts')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('K', p[0], p[1])
        return p[1]

    @_('BEGINNING duration', 'BEGINNING time',
       'BEGINNING person',)
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('BEGINNING', p[0], p[1])
        return p[1]

    @_('MIDDLE duration', 'MIDDLE time', 'MIDDLE num_dur')
    def simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('MIDDLE', p[0], p[1])
        return p[1]
 
    @_('Ø MIDDLE time')
    def simul(self, p):
        add_function('simultaneous', p[2])
        add_prep('MIDDLE', p[1], p[2])
        return p[2]
 
    @_('END duration', 'END time', 
       'END num_year', 'END num_day',
       'END person', 'END num_dur')
    def simul(self, p):
        conv_nums('NUMQ', p, as_dur=True)
        add_function('simultaneous', p[1])
        add_prep('END', p[0], p[1])
        return p[1]

    # -- CALENDRICAL TIMES -- 
    @_('B month')
    def month_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('B year', 'B num_year')
    def year_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('year_simul num_year')
    def year_simul(self, p):
        return init_ctime(
            p[0], p[1],
            functions=['simultaneous'],
        )

    @_('B day', 'B num_day')
    def day_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('B card')
    def day_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('Ø card')
    def day_simul(self, p):
        add_function('simultaneous', p[1])
        return p[1]

    @_('day_simul year_simul')
    def day_simul(self, p):
        add_qual('durative', p[1])
        return init_ctime(
            p[0], p[1],
            functions=['simultaneous'],
        )

    @_('B one_day')
    def oneday_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('B ordinal')
    def ordn_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('B first')
    def first_simul(self, p):
        add_function('simultaneous', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    # -- CALENDRICAL COMPOSITIONS --
    @_('day_simul month_simul', 
       'month_simul day_simul', 
       'month_simul year_simul',
       'month_simul year_simul day_simul',
       'month_simul year',
       'year_simul month_simul day_simul',
       'year_simul month_simul', 
        # NB: last month_simul below belongs to day_simul
        # this is a work-around, but might be a better model
       'year_simul month_simul day_simul month_simul',
       'year_simul month_simul oneday_simul',
       'year_simul day_simul',
    )
    def cal_simul(self, p):
        conv_nums('CALNUM', p)
        getattr(p, 'year', {})['functions'] =  ['simultaneous']
        return init_ctime(
            *p,
            functions=['cal_simul']
        )

    @_(
       'first_simul day_simul',
       'day_simul ordn_simul year_simul',
       'month_simul ordn_simul year_simul',
       'ordn_simul day_simul',
       'year_simul ordn_simul day_simul',
       'year_simul first_simul day_simul',
    )
    def cal_simul(self, p):
        try:
            ordtime = p.ordn_simul
        except AttributeError:
            ordtime = p.first_simul
        add_ref('CALORDN', ordtime['ordn'], ordtime)
        del ordtime['ordn']
        conv_nums('CALNUM', p)
        return init_ctime(
            *p,
            functions=['cal_simul'],
        )

    # -- either atelic ext or simul --
    @_('Ø one_time', 
       'Ø year', 'Ø month')
    def atelic_simul(self, p):
        add_function('atelic_ext, simultaneous', p[1])
        return p[1]
    
    @_('Ø one_day')
    def atelic_simul(self, p):
        add_function('simultaneous', p[1])
        conv_nums('NUMQ', p)
        return p[1]

    # -- telic extent / distances --
    @_('B num_dur', 'B one_time')
    def in_dur(self, p):
        conv_nums('NUMQ', p)
        add_function('telic_ext, dist_fut, dist_past', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    @_('B one_year')
    def in_dur(self, p):
        conv_num('NUMQ', p[1], as_dur=True)
        add_function('telic_ext, dist_fut, dist_past', p[1])
        add_prep('B', p[0], p[1])
        return p[1]

    # -- anterior --
    @_('L PNH/ duration', 
       'L PNH/ time', 
       'L PNH/ one_time',
       'L PNH/ all_dur',)
    def anterior(self, p):
        add_function('anterior', p[2])
        add_prep('PNH/', p[1], p[2])
        add_prep('L', p[0], p[2])
        return p[2]

    @_('L SFX PNH/')
    def anterior(self, p):
        time = init_time(p[2])
        add_function('anterior', time)
        add_ref('SFX', p[1], time)
        add_prep('L', p[0], time)
        return time

    @_('L PNH/ posts')
    def anterior(self, p):
        add_function('anterior', p[2])
        add_prep('L', p[0], p[2])
        add_prep('PNH/', p[1], p[2])
        return p[2]

    # stand-alone VRM/, will only match if 
    # followed by nothing else 
    @_('VRM/')
    def anterior(self, p):
        return init_time(
            p[0],
            functions=['anterior']
        )

    @_('atelic_ext anterior')
    def anterior_dist(self, p):
        return init_ctime(
            *p,
            functions=['anterior_dist']
        )

    # -- anterior durative --
    @_('<D duration', '<D time', 
       '<D simul', '<D one_time', 
       '<D year', '<D num_year', 
       '<D num_day', '<D one_day',
       '<D month', '<D person', '<D num_dur',
       '<D all_dur')
    def anterior_dur(self, p):
        conv_nums('NUMQ', p, as_dur=True)
        add_function('anterior_dur', p[1])
        add_prep('<D', p[0], p[1])
        add_qual('durative', p[1])
        return p[1]

    @_('L posterior')
    def anterior_dur(self, p):
        add_function('anterior_dur', p[1])
        add_prep('<D', p[0], p[1])
        add_qual('durative', p[1])
        return p[1]

    @_('L multi_time')
    def anterior_dur(self, p):
        add_function('anterior_dur', p[1])
        add_prep('L', p[0], p[1])
        add_qual('durative', p[1])
        p[1]['distributive'] = True
        return p[1]

    @_('<D multi_time')
    def anterior_dur(self, p):
        add_function('anterior_dur', p[1])
        add_prep('<D', p[0], p[1])
        add_qual('durative', p[1])
        p[1]['distributive'] = True
        return p[1]

    @_('anterior_dur anterior_dur')
    def multi_antdur(self, p):
        return init_ctime(
            *p,
            functions=['multi_antdur']
        )

    @_('antdur_simul anterior_dur')
    def multi_antdur(self, p):
        p[0]['functions'][0] = 'anterior_dur'
        add_qual('durative', p[0])
        return init_ctime(
            *p,
            functions=['multi_antdur']
        )

    # -- anterior dur / simultaneous
    @_('L time', 'L duration', 'L num_dur',
       'L num_day', 'L one_day',
       'L all_dur')
    def antdur_simul(self, p):
        conv_nums('NUMQ', p)
        add_function('anterior_dur, simultaneous', p[1])
        add_prep('L', p[0], p[1])
        return p[1]

    @_('L first')
    def antdur_simul(self, p):
        del p[1]['ordn']
        conv_nums('NUMQ', p)
        add_function('anterior_dur, simultaneous', p[1])
        add_prep('L', p[0], p[1])
        return p[1]

    @_('>L time')
    def antdur_simul(self, p):
        add_function('anterior_dur, simultaneous', p[1])
        add_prep('>L', p[0], p[1])
        return p[1]

    @_('L one_time')
    def simul(self, p):
        conv_num('NUMQ', p[1])
        add_function('simultaneous', p[1])
        add_prep('L', p[0], p[1])
        return p[1]

    @_('<D antdur_simul')
    def anterior_dur(self, p):
        p[1]['functions'][0] = 'anterior_dur'
        add_prep('<D', p[0], p[1])
        return p[1]

    @_('L posts', 'L posterior_dur')
    def antpost_dur(self, p):
        p[1]['functions'][0] = 'posterior_dur' 
        add_function('anterior_dur_past?', p[1])
        add_prep('L', p[0], p[1])
        return p[1] 

    @_('<D posts')
    def antpost_dur(self, p):
        p[1]['functions'][0] = 'posterior_dur' 
        add_function('anterior_dur_past?', p[1])
        add_prep('<D', p[0], p[1])
        return p[1] 

    @_('<D posterior')
    def antpost_dur(self, p):
        add_function('anterior_dur_past?', p[1])
        add_prep('<D', p[0], p[1])
        return p[1]

    # -- posteriors --
    @_('>XR/ duration', '>XR/ time', '>XR/ one_time',
       '>XR/ all_dur')
    def posterior(self, p):
        conv_nums('NUMQ', p)
        add_function('posterior', p[1])
        add_prep('>XR/', p[0], p[1])
        return p[1]

    @_('SFX >XR/')
    def posterior(self, p):
        time = init_time(p[1])
        add_function('posterior', time)
        add_ref('SFX', p[0], time)
        return time

    @_('>XR/ PERSON')
    def posterior(self, p):
        time = init_time(p[1])
        add_function('posterior', time)
        add_prep('>XR/', p[0], time) 
        add_ref('PERS', p[1], time)
        return time

    # stand-alone >XR; will only match if 
    # followed by nothing else 
    @_('>XR/')
    def posterior(self, p):
        return init_time(
            p[0], 
            functions=['posterior'],
        )

    @_('atelic_ext posterior')
    def posterior_dist(self, p):
        return init_ctime(
            p[0], p[1],
            functions= ['posterior_dist'],
        )

    # -- posteriors (durative?) --
    @_('MN time', 'MN one_time', 'MN one_day',
       'MN first', 'MN year', 'MN num_year', 'MN month',
       'MN day', 'MN dist_past', 'MN distago')
    def posts(self, p):
        conv_nums('NUMQ', p)
        add_function('posterior, posterior_dur', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    @_('MN duration', 'MN num_dur')
    def posts(self, p):
        add_function('posterior, posterior_dur', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    @_('MN multi_time')
    def posts(self, p):
        add_qual('durative', p[1])
        add_function('posterior, posterior_dur', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    @_('MN simul', 'MN first_simul')
    def posterior(self, p):
        add_function('posterior', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    @_('MN posterior')
    def posterior(self, p):
        add_function('posterior', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    # NB how the opposing duration leads
    # to a re-evaluation of this type!
    @_('MN anterior_dur')
    def posterior_dur(self, p):
        add_function('posterior_dur', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    @_('MN antdur_simul')
    def posterior_dur(self, p):
        add_function('posterior_dur', p[1])
        add_prep('MN', p[0], p[1])
        return p[1]

    # -- durations --
    @_('TIMES')
    def duration(self, p):
        return init_time(p[0], quals=['durative'])

    @_('DURATION')
    def duration(self, p):
        return init_time(p[0], quals=['durative'])

    @_('person_refm duration')
    def duration(self, p):
        return add_ref('PERS', p[0], p[1])

    @_('NUM TIMES', 'NUM_ONE TIMES')
    def num_dur(self, p):
        dur = init_time(
            p[1], quals=['durative'],
        )
        add_quant('NUMQ', p[0], dur)
        return dur 

    # ambiguity requires this tortuous approach; namely:
    # >> NUM TIMES NUM YEAR
    # if duration -> TIMES, and duration -> NUM duration
    # then the parser will shift duration without assigning
    # it to its number quantifier; Thus we spell out very
    # explicitly, which requires us to write separate patterns
    # for the reference tags
    @_('NUM THE TIMES', 'THE NUM TIME')
    def num_dur(self, p):
        dur = init_time(
            p[2], quals=['durative'],
        )
        add_ref('THE', p.THE, dur)
        add_quant('NUMQ', p.NUM, dur)
        return dur 

    @_('NUM MOD TIMES')
    def num_dur(self, p):
        dur = init_time(
            p[2], qualts=['durative'],
        )
        add_mod(p[1], dur)
        add_quant('NUMQ', p[0], dur)
        return dur

    @_('NUM GENREF TIMES')
    def num_dur(self, p):
        dur = init_time(
            p[2], quals=['durative'],
        )
        add_ref('GENREF', p[1], dur)
        add_quant('NUMQ', p[0], dur)
        return dur 

    @_('NUM time')
    def num_dur(self, p):
        add_quant('NUMQ', p[0], p[1])
        add_qual('durative', p[1])
        return p[1]

    @_('THIS num_dur')
    def num_dur(self, p):
        return add_ref('THIS', p[0], p[1])

    @_('ALL duration', 
       'ALL time',
       'ALL year')
    def all_dur(self, p):
        add_quant('ALL', p[0], p[1])        
        add_qual('durative', p[1])
        return p[1]

    @_('all_dur all_dur')
    def all_dur(self, p):
        if type(p[1]['times'][0]) == dict:
            return mergetime(
                p[0], p[1]
            )
        else:
            return init_ctime(
                *p,
                quals=['durative'], 
            )

    @_('ALL multi_time')
    def all_dur(self, p):
        add_qual('distributive', p[1])
        add_quant('ALL', p[0], p[1])
        return p[1]

    @_('MANY time', 'MANY duration')
    def duration(self, p):
        add_quant('MANY', p[0], p[1])        
        add_qual('durative', p[1])
        return p[1]

    @_('GENREF duration')
    def duration(self, p):
        return add_ref('GEN', p[0], p[1])

    @_('duration duration')
    def duration(self, p):
        return init_ctime(
            p[0], p[1],
            quals=['durative'],
        )

    @_('duration year',
       'duration month')
    def duration(self, p):
        time = init_ctime( 
            p[0], p[1],
            quals=['durative']
        )
        return time

    @_('num_dur num_year', 'num_day num_day',
       'num_dur num_day', 'num_year num_year',
       'num_day num_dur', 'num_dur num_dur')
    def num_dur(self, p):
        # set numbers to quantifiers
        conv_nums('NUMQ', p, as_dur=True)
        time = init_ctime( 
            p[0], p[1],
            quals=['durative']
        )
        return time

    @_('duration time', 'time duration')
    def duration(self, p):
        return init_ctime(
            p[0], p[1],
            quals=['durative'],
        )

    @_('time time', 'year year')
    def multi_time(self, p):
        return init_ctime(
            p[0], p[1],
        )

    @_('GENDUR time', 'GENDUR month', 'GENDUR duration')
    def duration(self, p):
        add_quant('DUR', p[0], p[1])
        add_qual('durative', p[1])
        return p[1]

    # add reference source_data
    @_('SFX duration')
    def duration(self, p):
        return add_ref('SFX', p[0], p.duration)

    @_('THE duration')
    def duration(self, p):
        return add_ref('THE', p[0], p.duration)

    @_('THIS duration')
    def duration(self, p):
        return add_ref('THIS', p[0], p.duration)

    @_('THAT duration')
    def duration(self, p):
        return add_ref('THAT', p[0], p.duration)
    
    # 2 Chr 21:15, ימים על־ימים
    # a kind of 'intensive' adjective
    @_('duration <L duration')
    def duration(self, p):
        add_prep('<L', p[1], p[2])        
        time = init_ctime(
            p[0], p[2],
            quals=['intensive'],
        )
        return time

    @_('year >XR/ year')
    def duration(self, p):
        add_prep('<L', p[1], p[2])        
        time = init_ctime(
            p[0], p[2],
            quals=['intensive'],
        )
        return time

    @_('NUM_ONE time')
    def one_time(self, p):
        add_quant('NUMQ', p[0], p[1])
        return p.time

    # -- time --
    @_('TIME')
    def time(self, p):
        return init_time(p.TIME)

    @_('TIME2')
    def time(self, p):
        return init_time(p.TIME2, duplicate=True)
   
    @_('SFX time')
    def time(self, p):
        return add_ref('SFX', p[0], p.time)

    @_('THE time')
    def time(self, p):
        return add_ref('THE', p[0], p.time)

    @_('THIS time')
    def time(self, p):
        return add_ref('THIS', p[0], p.time)

    @_('THAT time')
    def time(self, p):
        return add_ref('THAT', p[0], p.time)

    @_('ORDNM time')
    def time(self, p):
        return add_ref('ORDN', p[0], p[1])

    @_('person_refm time')
    def time(self, p):
        return add_ref('PERS', p[0], p[1])

    @_('MOD time')
    def time(self, p):
        return add_mod(p[0], p[1])

    @_('GENREF time')
    def time(self, p):
        return add_ref('GEN', p[0], p[1])

    # define specialized time
    @_('BEGINNING', 'END')
    def time(self, p):
        return init_time(p[0])

    # calendrical time references
    @_('DAY')
    def day(self, p):
        return init_time(p.DAY)

    @_('THE day')
    def day(self, p):
        return add_ref('THE', p[0], p.day)

    @_('ORDNM day')
    def day(self, p):
        return add_ref('ORDN', p[0], p[1])

    @_('GENREF day')
    def day(self, p):
        return add_ref('GEN', p[0], p[1])

    @_('NUM_ONE day')
    def one_day(self, p):
        p.day['NUM'] = p[0]
        check_nums(p.day)
        return p.day

    @_('NUM day')
    def num_day(self, p):
        p.day['NUM'] = p[0]
        check_nums(p.day)
        return p.day

    @_('THIS num_day')
    def num_day(self, p):
        return add_ref('THIS', p[0], p[1])

    @_('num_day month_ref', 
       'num_day person_ref',
       'one_day month_ref')
    def day(self, p):
        conv_num('CALNUM', p[0])   
        return add_ref(
            p[1]['ref'], 
            p[1],
            p[0],
        )

    @_('day month_ref')
    def day(self, p):
        return add_ref(
            p[1]['ref'],
            p[1],
            p[0],
        )

    @_('CARD')
    def card(self, p):
        return init_time(
            p[0],
            NUM=p[0],
        )

    @_('CARDC card')
    def card(self, p):
        # replace single number
        # with whole card chain
        # this is a workaround since 
        # the "head" of the cardinal chain
        # is selected by default
        return init_time(p[0]) 

    @_('THE card')
    def card(self, p):
        add_ref('THE', p[0], p[1])
        return p[1]

    @_('card month_ref', 'card ordn_ref', 
       'day ordn_ref')
    def day(self, p):
        conv_nums('CALNUM', p)
        ref = p[1]['ref']
        del p[1]['ref']
        if 'ordn' in p[1]:
            del p[1]['ordn']
        return add_ref(
            ref,
            p[1],
            p[0],
        )

    @_('MONTH')
    def month(self, p):
        return init_time(p[0])

    @_('THE month')
    def month(self, p):
        return add_ref('THE', p[0], p[1])

    @_('THIS month')
    def month(self, p):
        return add_ref('THIS', p[0], p[1])

    @_('THAT month')
    def month(self, p):
        return add_ref('THAT', p[0], p[1])

    @_('ORDNM month')
    def month(self, p):
        return add_ref('ORDN', p[0], p[1])

    @_('GENREF month')
    def month(self, p):
        return add_ref('GEN', p[0], p[1])

    @_('NUM month')
    def month(self, p):
        return add_ref('CALNUM', p[0], p[1])

    @_('SFX month')
    def month(self, p):
        return add_ref('SFX', p[0], p[1])

    @_('month antdur_simul')
    def month(self, p):
        p[1]['functions'][0] = 'reference'
        return add_ref('TIME', p[1], p[0]) 

    @_('month L year')
    def month(self, p):
        add_function('reference', p[2])
        add_prep('L', p[1], p[2])
        return add_ref('YEAR', p[2], p[0])

    @_('GENCARD day')
    def day(self, p):
        return add_ref('CALNUM', p[0], p[1])

    @_('YEAR')
    def year(self, p):
        return init_time(p.YEAR)

    @_('THE year')
    def year(self, p):
        return add_ref('THE', p[0], p[1])

    @_('THIS year') 
    def year(self, p):
        return add_ref('THIS', p[0], p[1])

    @_('THAT year')
    def year(self, p):
        return add_ref('THAT', p[0], p[1])

    @_('ORDNM year')
    def year(self, p):
        return add_ref('ORDN', p[0], p[1])

    @_('MOD year')
    def year(self, p):
        return add_mod(p[0], p[1])

    @_('GENREF year')
    def year(self, p):
        return add_ref('GEN', p[0], p[1])

    @_('GENCARD year')
    def year(self, p):
        return add_ref('CALNUM', p[0], p[1])

    @_('NUM year')
    def num_year(self, p):
        p.year['NUM'] = p.NUM
        check_nums(p[1])
        return p.year

    @_('THIS num_year')
    def num_year(self, p):
        return add_ref('THIS', p[0], p[1])

    @_('NUM_ONE year')
    def one_year(self, p):
        p.year['NUM'] = p.NUM_ONE
        check_nums(p[1])
        return p.year

    @_('year person_ref')
    def year(self, p):
        return add_ref('PERS', p[1], p.year) 
    
    @_('num_year person_ref')
    def year(self, p):
        conv_num('CALNUM', p.num_year)
        add_ref('PERS', p[1], p.num_year)
        return p.num_year
    
    @_('year antdur_simul')
    def year(self, p):
        p[1]['functions'][0] = 'reference'
        return add_ref('TIME', p[1], p[0]) 

    # -- adverb time --
    @_('NOW')
    def time(self, p):
        time = init_time(p[0])
        add_tense('PRES', p[0], time)
        add_ref('DEICTIC', p[0], time)
        return time

    @_('TOMORROW')
    def time(self, p):
        time = init_time(p[0])
        add_tense('FUT', p[0], time)
        add_ref('DEICTIC', p[0], time)
        return time

    @_('YESTERDAY')
    def time(self, p):
        time = init_time(p[0])
        add_tense('PAST', p[0], time)
        add_ref('DEICTIC', p[0], time)
        return time

    @_('DISTAGO')
    def distago(self, p):
        time = init_time(p[0])
        add_tense('PAST', p[0], time)
        add_quant('DISTAGO', p[0], time)
        add_ref('DEICTIC', p[0], time)
        return time

    @_('THUS')
    def time(self, p):
        time = init_time(p[0])
        return add_ref('DEICTIC', p[0], time)

    @_('THEN')
    def time(self, p):
        time = init_time(p[0])
        return add_ref('DEICTIC', p[0], time)

    @_('ORDN')
    def ordinal(self, p):
        time = init_time(p[0])
        add_ref('ORDN', p[0], time)
        time['ordn'] = p.ORDN
        return time

    @_('THE ordinal')
    def ordinal(self, p):
        return add_ref('THE', p[0], p[1]) 

    # see Gen 28:19 for why this pattern is needed:
    # לראשנה
    # to treat RISHON as a normal ordinal would mean
    # that this phrase would get interpreted as a 
    # calendrical reference; whereas it never occurs
    # as such in calendrical simultaneous phrases
    @_('FIRST')
    def first(self, p):
        time = init_time(p[0])
        time['ordn'] = p[0]
        return time

    @_('THE first')
    def first(self, p):
        return add_ref('THE', p[0], p[1])

    @_('PERSON')
    def person(self, p):
        pers = init_time(p[0]) 
        add_ref('PERS', p[0], pers)
        return pers

    @_('THE person')
    def person(self, p):
        add_ref('THE', p[0], p[1])
        return p[1]

    @_('SFX person')
    def person(self, p):
        add_ref('SFX', p[0], p[1])
        return p[1]

    @_('GENREF person')
    def person(self, p):
        return add_ref('GEN', p[0], p[1])
    
    # NB: reference constructions are not 
    # treated as typical times; so we convert
    # the time to a simplified dictionary
    @_('L person')
    def person_ref(self, p):
        add_prep('L', p[0], p[1])
        add_function('reference', p[1])
        p[1]['ref'] = 'PERS'
        return p[1]

    @_('PERSONM')
    def person_refm(self, p):
        return init_time(
            p[0],
            functions=['reference'],
            ref='PERS'
        )

    @_('L month')
    def month_ref(self, p):
        add_prep('L', p[0], p[1])
        add_function('reference', p[1])
        p[1]['ref'] = 'MONTH'
        return p[1]

    @_('MONTHREF month')
    def month_ref(self, p):
        add_function('reference', p[1])
        p[1]['ref'] = 'MONTH'
        return p[1]

    @_('IT B')
    def month_ref(self, p):
        time = init_time(p[1])
        add_ref('SFX', p[0], time)
        add_prep('B', p[1], time)
        add_function('simultaneous', time)
        add_function('reference', time)
        time['ref'] = 'MONTH'
        return time 

    @_('L ordinal')
    def ordn_ref(self, p):
        add_prep('L', p[0], p[1])
        add_function('reference', p[1])
        p[1]['ref'] = 'CALORDN'
        return p[1]

def parse_times(paths, API):
    """Apply the time parser."""

    # load phrase parsings
    phrases = ParseLoader(paths['ph_parses']).load()
    functs = ParseLoader(paths['functions']).load()

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
    
    # gather eligible clauses;
    # an eligible clause is one for which all
    # its Time phrases are parsed; There is some
    # complexity here in that some phrases have been
    # subsumed into another's parse. We use the 
    # value "has_phrases" to inventory all phrases that are
    # parsed; then we cross reference that with the time
    # phrases found in the clause; finally we filter out 
    # all of the clause's time phrases that have indeed
    # been subsumed since they will be parsed as part of
    # their head phrase 
    clauses = {}
    parsed_phrases = set(
        ph_node for head_ph in phrases
            for ph_node in phrases[head_ph]['has_phrases']
    )
    for clause in F.otype.s('clause'):
        time_phrases = [
            ph for ph in L.d(clause, 'phrase')
                if functs[ph] == 'Time'
        ]
        if time_phrases and parsed_phrases.issuperset(set(time_phrases)):
            # filter out those phrases which are subsumed
            clauses[clause] = [
                ph for ph in time_phrases 
                    if ph in phrases
            ]

    print(len(clauses), 'clauses ready for parsing...')

    for clause, time_phrases in clauses.items(): 

#        if clause == 427993:
#            tokenizer.debug = True
#        else:
#            tokenizer.debug = False

        # collect tokens for all phrase parsings
        tokens = []
        for ph in time_phrases:
            parsing = phrases[ph]['parse']
            ph_tokens = tokenizer.tokenize(parsing)
            tokens.extend(
                ph_tokens
            )

        # attempt parsing
        parsing = parser.parse(t for t in tokens)
        if parsing is not None and error_tracker['e'] is None:
            parsing['phrase_nodes'] = time_phrases
            parsed[clause] = parsing
        else:
            toks = str([t.type for t in tokens])
            errors[clause] = (error_tracker['e'], toks)

            # reset error tracker
            error_tracker['e'] = None 

    # export
    with open(paths['parsed'], 'w') as outfile:
        json.dump(
            parsed, outfile, 
            indent=2, ensure_ascii=False,
        )

    with open(paths['notparsed'], 'w') as outfile:
        json.dump(errors, outfile, indent=2, ensure_ascii=False) 
