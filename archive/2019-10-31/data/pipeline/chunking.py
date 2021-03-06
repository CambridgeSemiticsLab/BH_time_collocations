import collections
from tf.fabric import Fabric

class Chunker:
    '''
    Class performs three main actions:
        1. constructs timephrase chunks by merging
           split adjacent time phrases + split modifiers
        2. chunks quantifier chains
        3. chunks preposition chains
    These objects are built on top of BHSA objects w/
    modified phrase function data.
    All data is exported to the export directory.
    '''
    
    def __init__(self, locs, metadata):
        '''
        At initialization, build node counter that calculates maxNode
        from BHSA. This counter is operated on by all methods. Init also
        builds the dictionaries where all new node and edge data is 
        stored.
        '''
        
        # NB that new function data is loaded from output_dir
        self.output_dir = locs['output']
        locas = [locs['bhsa'], locs['heads'], locs['output']]
        TF = Fabric(locations=locas, silent=True)
        load_features = ['function', 'note', 'sem_set', 'head', 'nhead',
                         'st', 'ls', 'language', 'pdp', 'lex', 'obj_prep', 'language']
        self.bhsa = bhsa = TF.load(' '.join(load_features), silent=True)
        self.meta = metadata
        
        # configure counters and dicts to build/store new objects
        self.nodeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        self.edgeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        
        # gather BHSA otype and oslots
        # new objects/features will be generated on top of these
        self.nodeFeatures['otype'] = dict((n, bhsa.F.otype.v(n)) for n in bhsa.N())
        self.edgeFeatures['oslots'] = dict((n, bhsa.L.d(n, 'word')) for n in bhsa.N() 
                                               if bhsa.F.otype.v(n) != 'word')
        
        # new node numbers calculated from here
        self.newNode = max(self.nodeFeatures['otype'].keys())
        
    def chunk_time(self):
        '''
        There are several cases in BHSA where time phrases 
        are divided into several pieces whereas elsewhere 
        the parts are kept together as a single phrase. 
        This is an undesirable inconsistency. To solve this problem, 
        new chunk objects are generated and mapped over the two phrase
        boundaries. The objects have a label=timephrase.
        This method operates on self.edgeFeatures and self.nodeFeatures
        to add new objects.
        '''
        
        # define shortform TF methods for easy use
        F, E, L = self.bhsa.F, self.bhsa.E, self.bhsa.L
        
        # iterate through all phrases in BHSA
        # if phrase function=Time and is not followed by 
        # another time phrase, generate a timephrase chunk.
        # Otherwise, if it is followed, include the subsequent
        # time phrase in the chunk. Skip over times that are preceded
        # by another time since these will be subsumed into the first time.
        for phrase in F.function.s('Time'):
            
            # only chunk Hebrew cases
            language = F.language.v(L.d(phrase, 'word')[0])
            if language != 'Hebrew':
                continue
            
            # assign chunk boundaries here
            chunkSlots = []
            
            # examine positions around time phrase
            thisclause = L.d(L.u(phrase, 'clause')[0], 'phrase') # all phrases in cl
            nextphrase = next(iter(L.n(phrase, 'phrase')), 0)
            prevphrase = next(iter(L.p(phrase, 'phrase')), 0)
            nx_time = F.function.v(nextphrase) == 'Time' and nextphrase in thisclause
            pr_time = F.function.v(prevphrase) == 'Time' and prevphrase in thisclause
            
            # Automatically chunk time phrases not preceded 
            # or followed by another time prhase
            if not nx_time and not pr_time:
                chunkSlots = L.d(phrase, 'word')
            
            # Chunk all first position times and their
            # subsequent daughters into a single chunk
            if nx_time and not pr_time:
                
                chunkSlots.extend(L.d(phrase, 'word')) # count existing slots

                # gather all subsequent slots
                # iteratively reassign nx_time to the next phrase
                # until all daughters are captured
                while F.function.v(nextphrase) == 'Time' and nextphrase in thisclause:
                    chunkSlots.extend(L.d(nextphrase, 'word'))
                    nextphrase = next(iter(L.n(nextphrase, 'phrase')), 0)
                
            # skip non-dominant time phrases
            elif pr_time:
                continue
            
            # finalize time chunk object
            self.newNode += 1
            self.edgeFeatures['oslots'][self.newNode] = chunkSlots
            self.nodeFeatures['otype'][self.newNode] = 'timephrase'
            # add semantic roles for time heads
            for time in E.nhead.t(phrase):
                self.edgeFeatures['role'][time] = {self.newNode:'time'}
    
    def get_item(self, item, iterable, default=0):
        '''
        Checks if item is in an iterable.
        If present, return item. If absent, return default.
        Good for selecting +/- slots in a clause/phrase tuple.
        '''
        if item in iterable:
            return item
        else:
            return default
        
    def fill_slots(self, chunk):
        '''
        Fills in gapped slots such as waws and other
        items that are missing in a quant chunk (see below).
        '''
        chunk.sort()
        minSlot, maxSlot = chunk[0], chunk[-1]
        return list(range(minSlot, maxSlot+1))
        
    def chunk_quants(self):
        '''
        In Hebrew, cardinal numbers are split
        into several parts, e.g. חמשׁ ועשרים יום
        == "five and twenty days" For this project, 
        I want to be able to easily select and label 
        number chains. This method does that by 
        chunking numbers into chunk objects.         
        The chunks are built in three steps:
        1. quantifier "atoms" are isolated, which consist of
            the following pairs:
            * cardinal + noun
            * noun + cardinal
            * cardinal + ה + noun
            * cardinal + ו + cardinal
            * cardinal + cardinal
            These pairs account for all of the atomic parts
            of a quantifier chunk.
        2. combine quantifier atoms by isolating continuous
            runs. These are basic quant chunks.
        3. combine adjacent quant chunks to capture larger,
            composite chains: e.g. שבע שנה ומאת שנה.
        '''
        
        # define shortform TF methods for easy use
        F, E, L = self.bhsa.F, self.bhsa.E, self.bhsa.L
        
        # Gather qantifier atoms
        quant_atoms = []
        cardinals = [word for lex in F.ls.s('card')
                        for word in L.d(lex, 'word')
                        if F.language.v(word) == 'Hebrew']
        
        for number in cardinals:
            
            # get positions above and around number
            phrase = L.d(L.u(number, 'phrase')[0], 'word')
            phrase_atom = L.d(L.u(number, 'phrase_atom')[0], 'word')
            p1 = self.get_item(number+1, phrase, 0) # plus 1 slot, etc.
            p2 = self.get_item(number+2, phrase, 0) 
            m1 = self.get_item(number-1, phrase, 0) # minus 1, etc.
            m2 = self.get_item(number-2, phrase, 0)
            m1pa = self.get_item(number-1, phrase_atom, 0) # minus 1 in phrase atom
            subs = {'subs', 'nmpr', 'prps'} # valid POS for quantified objects
            
            # match atoms
            # card + noun
            if all([F.pdp.v(p1) in subs,
                    F.ls.v(p1) != 'card',
                    F.sem_set.v(p1) != 'prep']):
                quant_atoms.append((number, p1))
                
            # noun + card
            elif all([F.pdp.v(m1pa) in subs,
                      F.ls.v(m1pa) != 'card',
                      F.sem_set.v(m1pa) != 'prep',
                      F.st.v(m1pa) == 'a']):
                quant_atoms.append((m1pa, number))
            
            # card + ה + noun
            elif all([F.lex.v(p1) == 'H',
                      F.pdp.v(p2) in subs,
                      F.ls.v(p2) != 'card']):
                quant_atoms.append((number, p2))
                
            # card + ו + card
            elif all([F.lex.v(p1) == 'W',
                      F.ls.v(p2) == 'card']):
                quant_atoms.append((number, p2))
                
            # card + card
            elif F.ls.v(p1) == 'card':
                quant_atoms.append((number, p1))
                
            # ø + card + ø
            elif not any([F.lex.v(m1) == 'W' and F.ls.v(m2) == 'card', 
                          F.ls.v(m1) == 'card', 
                          F.lex.v(number)=='>XD/']):
                quant_atoms.append((number,))
                
        quant_atoms = sorted(quant_atoms) # sort chunks for processing
                
        # link atoms on the last number in the chains
        # e.g. (3800, 3802), (3802, 3803), (3803, 3804)
        # links should be made with 3802 and 3803, which
        # unite the three tuples. We do this with a dict.
        # The key is the last item in the tuple, the value
        # is a list of the tuple's contents. 
        # We use dict.pop() when a tuple's first
        # item matches another's last item, the key, 
        # and we pop the value of the first & combine the 2
        qchunks = {}
        for atom in quant_atoms:
            # match first slot to last slot in an atome
            if atom[0] in qchunks:
                qchunks[atom[-1]] = qchunks.pop(atom[0]) + list(atom)[1:]
            # map atoms to dict
            else:
                qchunks[atom[-1]] = list(atom)
        
        qchunks = list(qchunks.values()) # take the chunks
        
        # link chunks contained in same phrase_atom
        # this includes cases where two nouns are quantified
        # but both quant chunks are intended to be read together
        # as a sigle set. Often the two nouns are the same, but not
        # always, e.g. "forty days and forty nights."
        
        # make phrase_atom 2 chunk mapping
        phrase2chunks = collections.defaultdict(list)
        for chunk in qchunks:
            phrase_atom = L.u(chunk[0], 'phrase_atom')[0]
            phrase2chunks[phrase_atom].append(chunk)
            
        # find composite chunks
        for phrase, chunks in phrase2chunks.items():
            # skip simple chunks
            if len(chunks) < 2:
                continue
            # get the quantified nouns
            chunknouns = [word for chunk in chunks for word in chunk
                              if F.ls.v(word) != 'card']
            # look for composite quantification
            if len(chunknouns) < 2:
                continue
            # add composite chunk
            qchunks.append([word for chunk in chunks for word in chunk])
        
        # generate the chunk objects
        for chunk in qchunks:
            chunkslots = self.fill_slots(chunk) # fill in gaps
            chunknouns = [word for word in chunk if F.ls.v(word) != 'card']
            chunkquants = [word for word in chunk if F.ls.v(word) == 'card']
            label = 'quant_NP' if chunknouns else 'quant' # 2 labels depending on presence of noun(s)
            self.newNode += 1
            self.edgeFeatures['oslots'][self.newNode] = chunkslots
            self.nodeFeatures['otype'][self.newNode] = 'chunk'
            self.nodeFeatures['label'][self.newNode] = label
            self.edgeFeatures['role'].update({noun:{self.newNode:'subs'} for noun in chunknouns})
            self.edgeFeatures['role'].update({quant:{self.newNode:'quant'} for quant in chunkquants})
    
    def climb_prep_chain(self, prep, prep_list):
        '''
        Recursively climbs a prepositional chain (see next).
        '''
        # define shortform TF methods for easy use
        F, E = self.bhsa.F, self.bhsa.E
        prep_list.append(prep)
        daughter = next((po for po in E.obj_prep.t(prep) if F.sem_set.v(po)=='prep'),[])
        if daughter:
            self.climb_prep_chain(daughter, prep_list)
    
    def chunk_preps(self):
        '''
        Prepositions that are chained together function 
        as a single directional unit, and some words function 
        as prepositions within a certain frame where elsewhere 
        they may function as nouns. Using the sem_set feature 
        from the heads project and the obj_prep edge relation, 
        we can easily export a construction that can cover these cases.
        '''
        
        # define shortform TF methods for easy use
        F, E = self.bhsa.F, self.bhsa.E
        
        for prep in F.sem_set.s('prep'):
            
            # skip subordinated preps
            if E.obj_prep.f(prep):
                continue
                
            # climb down prep chain
            prep_cx = []
            self.climb_prep_chain(prep, prep_cx)

            # export object
            self.newNode += 1
            self.edgeFeatures['oslots'][self.newNode] = prep_cx
            self.nodeFeatures['otype'][self.newNode] = 'chunk'
            self.nodeFeatures['label'][self.newNode] = 'prep'
    
    def execute(self):
        '''
        Runs the chunking methods in the 
        necessary succession and exports the
        new TF data.
        '''
        
        print('Running chunkers...')
        print('\trunning time chunker...')
        self.chunk_time()
        print('\trunning quant chunker...')
        self.chunk_quants()
        print('\trunning prep chunker...')
        self.chunk_preps()
        
        # get report on all new objects
        chunks = [node for node, otype in self.nodeFeatures['otype'].items()
                    if otype=='chunk']
        
        nlabels = collections.Counter(self.nodeFeatures['label'][node] 
                                          for node in chunks)
        nroles = collections.Counter(role for node, edge in self.edgeFeatures['role'].items()
                                        for role in edge.values())

        print(f'{len(chunks)} chunk objects formed...')
        for label, count in nlabels.most_common():
            print(f'\t{count} chunk objects with label {label}')
        for label, count in nroles.most_common():
            print(f'\t{count} words with chunk role of {label}')
        
        
        # update metadata and export
        
        metadata = {'':self.meta,
                    'oslots': {'valueType':'int', 
                               'edgeValues':False},
                    'otype':{'valueType':'str'},
                    'role':{'edgeValues':True,
                            'valueType':'str', 
                            'description':'role of a word in a chunk object'},
                    'label':{'valueType':'str',
                             'description':'a label for a chunk object'},
                   }
        print('exporting tf...')
        TFexport = Fabric(locations=self.output_dir, silent=True)
        TFexport.save(metaData=metadata, nodeFeatures=self.nodeFeatures, edgeFeatures=self.edgeFeatures)
        print('SUCCESS')
