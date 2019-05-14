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
    
    def __init__(self, bhsa_dir, heads_dir, output_dir, metadata):
        '''
        At initialization, build node counter that calculates maxNode
        from BHSA. This counter is operated on by all methods. Init also
        builds the dictionaries where all new node and edge data is 
        stored.
        '''
        
        # NB that new function data is loaded from output_dir
        locas = [bhsa_dir, heads_dir, output_dir, metadata]
        TF = Fabric(locations=locas, silent=True)
        load_features = ['function', 'note', 'sem_set', 'head', 'nhead']
        self.bhsa = bhsa = self.TF.load(' '.join(load_features))        
        
        # configure counters and dicts to build/store new objects
        self.nodeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        self.edgeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        
        # build BHSA otype, oslots, and new function & note features
        # new objects will be generated on top of these
        # but all are needed in order to export the new oslots.tf and otype.tf files
        self.nodeFeatures['otype'] = dict((n, bhsa.F.otype.v(n)) for n in bhsa.N())
        self.edgeFeatures['oslots'] = dict((n, bhsa.L.d(n, 'word')) for n in bhsa.N()\ 
                                               if bhsa.F.otype.v(n) != 'word')
        self.nodeFeatures['function'] = dict((n, bhsa.F.function.v(n))\
                                                 for n in bhsa.F.otype.s('phrase'))
        self.nodeFeatures['note'] = dict((n, bhsa.F.note.v(n))\
                                             for n in bhsa.N() if bhsa.F.note.v(n))
        
        # new node numbers calculated from here
        self.newNode = max(self.nodeFeatures['otype'].keys())
        
        
    def chunk_time(self):
        '''
        There are several cases in BHSA where time phrases 
        are divided into several pieces whereas elsewhere 
        the parts are kept together as a single phrase. 
        This is an undesirable inconsistency. To solve this problem, 
        new phrase boundaries are generated and mapped over the old 
        boundaries stored in the oslots file. 
        
        This method operates on self.edgeFeatures and self.nodeFeatures
        to add new objects.
        '''
        
        # define shortform TF methods for easy use
        F, L = self.bhsa.F, self.bhsa.L
        
        # iterate through all phrases in BHSA
        # if phrase function=Time and is not followed by 
        # another time phrase, generate a timephrase chunk.
        # Otherwise, if it is followed, include the subsequent
        # time phrase in the chunk. Skip over times that are preceded
        # by another time since these will be subsumed into the first time.
        for phrase in F.function.s('Time'):
            
            # assign features if needed
            chunkSlots = []
            chunkRoles = {} # wordnode that is time mapped to this chunk
            
            # examine positions around time phrase
            thisclause = L.d(L.u(phrase, 'clause')[0], 'phrase') # all phrases in cl
            nextphrase = next(iter(L.n(phrase, 'phrase')), 0)
            prevphrase = next(iter(L.p(phrase, 'phrase')), 0)
            nx_time = F.function.v(nextphrase) == 'Time' and nx_time in thisclause
            pr_time = F.function.v(prevphrase) == 'Time' and pr_time in thisclause
            
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
                while F.function.v(nx_time) == 'Time' and nx_time in thisclause:
                    chunkSlots.extend(L.d(nx_time, 'word'))
                    nx_time = next(iter(L.n(nx_time, 'phrase')), 0)
                
            # skip non-dominant time phrases
            elif pr_time:
                continue
            
            # finalize time chunk object
            self.newNode += 1
            self.edgeFeatures['oslots'][self.newNode] = chunkSlots
            self.nodeFeatures['otype'][self.newNode] = 'chunk'
            self.nodeFeatures['label'][self.newNode] = 'timephrase'
            # add semantic roles for time heads
            for time in E.head.t(phrase):
                self.edgeFeatures['role'][time] = {phrase:'time'}
    
        
    def chunk_quants(self):
        '''
        In Hebrew, cardinal numbers are split
        into several parts, e.g. 
        חמשׁ ועשרים יום
        == "five and twenty days"
        For this project, I want to be able to easily
        navigate and filter through number chains.
        This module does that by chunking numbers into
        individual objects. That allows more complex chunks
        to be handled easily.
        '''
        
    
        
    def execute(self):
        '''
        Orders and runs all chunking methods
        '''
        
    