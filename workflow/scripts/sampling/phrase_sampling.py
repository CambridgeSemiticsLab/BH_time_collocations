import collections
import json
from tf.fabric import Fabric

def get_phrase_samples(tf_locs, sample_out, unsample_out):
    """Generate phrase samples using project directives."""
    
    # initialize TF with supplied paths
    TF = Fabric(locations=tf_locs)
    API = TF.load('language typ st function pdp')
    F, L = API.F, API.L

    # loop through all phrases in BHSA and collect samples
    samples = {}
    nonsamples = collections.defaultdict(list)
    for phrase_atom in F.otype.s('phrase_atom'):

        # skip ineligible phrase types
        keep_typs = {
            'NP', 'PP', 'AdjP', 'AdvP', 
            'PPrP', 'PrNP'
        }
        if F.typ.v(phrase_atom) not in keep_typs:
            nonsamples['typ'].append(phrase_atom)
            continue

        # skip ineligible functions
        phrase = L.u(phrase_atom, 'phrase')[0]
        keep_functs = {
            'Adju', 'Cmpl', 'Loca', 
            'Modi', 'Objc', 'PreC',
            'Subj', 'Time'
        }
        if F.function.v(phrase) not in keep_functs:
            nonsamples['function'].append(phrase_atom)
            continue

        # look up contextual data for filtering
        phrase_words = L.d(phrase, 'word')
        words = L.d(phrase_atom, 'word')

        # skip non-Hebrew phrase_atoms
        if F.language.v(words[0]) != 'Hebrew':
            nonsamples['aramaic'].append(phrase_atom)
            continue

        # Skip phrase_atoms contained in phrases that 
        # have interruptions from other embeddings.
        # consecutive slots:
        #     s0 + len(slots)-1 == s-1
        # where s0 is first slot in set and s-1 is last slot
        consec_end = phrase_words[0] + len(phrase_words)-1
        if consec_end != phrase_words[-1]:
            nonsamples['split'].append(phrase_atom)
            continue

        # skip phrase_atoms which have a word ending in construct
        if F.st.v(words[-1]) == 'c':
            nonsamples['C$'].append(phrase_atom)
            continue

        # anything that makes it here is a match
        samples[phrase_atom] = words

    # export the samples and non-samples
    with open(sample_out, 'w') as outfile:
        json.dump(samples, outfile, indent=1)

    with open(unsample_out, 'w') as outfile:
        json.dump(nonsamples, outfile, indent=1)
