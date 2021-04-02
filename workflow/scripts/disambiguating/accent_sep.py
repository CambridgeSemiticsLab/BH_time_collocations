"""
Find substantives that should have a separator 
between them and a subsequent substantive according
to the Tiberian accents.

One can think of it like an "invisible" conjunction.
These patterns are especially frequent in lists.
For example:

    "Shem Ham and Japheth"

Here we want to insert a separator between "Shem" and
"Ham" so that the parser knows these are separate and
not an apposition relation. This requires semantic knowledge.

We harvest some of this knowledge from the Tiberian accents.
First we find all such cases of side-by-side substantives.
Then we identify whether the first substantive has a disjunctive
accent or not. If it does, we add it to a list.
"""

import json
from .accents import AccentTagger

subs_query = '''
phrase_atom
    word pdp=nmpr|subs ls#card st=a
    <: word pdp=nmpr|subs ls#card
'''

def get_accent_seps(paths, tf_api):
    """Harvest accent separators from Tiberian accents."""

    # tags accents on BHSA words
    tagger = AccentTagger(tf_api)

    # execute the query to return BHSA nodes
    query = tf_api.S.search(subs_query)

    # gather words with disjunctive accents
    accent_seps = []
    for ph, w1, w2 in sorted(query):
        accent_type = tagger.tag(w1)
        if accent_type == 'disjunct':
            accent_seps.append(w1)

    # export results
    with open(paths['accentseps'], 'w') as outfile:
        json.dump(accent_seps, outfile)
