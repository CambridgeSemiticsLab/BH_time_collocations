import collections

def chunk_time(addTypes, api):
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
    F, E, L = api.F, api.E, api.L

    # iterate through all phrases in BHSA
    # if phrase function=Time and is not followed by 
    # another time phrase, generate a timephrase chunk.
    # Otherwise, if it is followed, include the subsequent
    # time phrase in the chunk. Skip over times that are preceded
    # by another time since these will be subsumed into the first time.
    oid = 1
    addTypes['timephrase'] = {
        'nodeSlots': collections.defaultdict(set),
    }
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
        addTypes['timephrase']['nodeSlots'][oid] = set(chunkSlots)
        oid += 1
