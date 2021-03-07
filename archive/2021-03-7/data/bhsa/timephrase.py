import collections
from positions import Walker

def chunk_time(addTypes, tf):
    """Combine cleft BHSA time phrases into a single unit. 

    There are several cases in BHSA where adjacent 
    time phrases are divided into several pieces 
    whereas elsewhere the parts are kept together 
    as a single phrase. This is an undesirable inconsistency. 
    To solve this problem, new chunk objects are generated 
    and mapped over the two phrase boundaries. The objects 
    have a label=timephrase. This method operates on self.edgeFeatures 
    and self.nodeFeatures to add new objects.

    Args:
        addTypes: a dictionary of otype to slots and feature mappings
            in the format of 
            >> addTypes[otype] = {nodeSlots: {node:slots_tuple}}
        tf: an instance of Text-Fabric with BHSA
    """

    # define shortform TF methods for easy use
    F, E, L = tf.api.F, tf.api.E, tf.api.L

    # validator function for identifying times
    def is_time(node):
        return F.function2.v(node) == 'Time'

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
    covered = set() # track already-processed phrases
    for phrase in F.function2.s('Time'):

        # skip covered phrases
        if phrase in covered:
            continue

        # only chunk Hebrew cases
        language = F.language.v(L.d(phrase, 'word')[0])
        if language != 'Hebrew':
            continue

        # assign chunk boundaries here
        chunkSlots = []

        # collect phrases within clause
        # tuple of phrases will serve as a path to walk
        # NB: we use a method of node-adjacency rather than Text-Fabric 
        # slot adjacency to capture gapped time phrases, 
        # the material in between will need to be captured later
        clause_phrases = L.d(L.u(phrase, 'clause')[0], 'phrase') # all phrases in cl
        Wk = Walker(phrase, clause_phrases) # walks clause phrases when called

        # collect all adjacent times in clause
        # Walker takes a function to validate nodes along the path
        # another function, stop, ends the walk 
        # results are returned as long as they are validated by val function
        times = [phrase]
        times.extend(
            Wk.ahead(
                is_time, 
                every=True, 
                stop=lambda n: not is_time(n),
                default=[],
            )
        )
        covered |= set(times)

        # build timephrase oslots
        for time in times:
            chunkSlots.extend(L.d(time,'word'))

        # finalize time chunk object
        addTypes['timephrase']['nodeSlots'][oid] = set(chunkSlots)
        oid += 1
