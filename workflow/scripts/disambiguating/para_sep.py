"""
Find substantives that should have a separator 
between them and a subsequent substantive according
to the ETCBC's BHSA.

One can think of it like an "invisible" conjunction.
These patterns are especially frequent in lists.
For example:

    "Shem Ham and Japheth"

Here we want to insert a separator between "Shem" and
"Ham" so that the parser knows these are separate and
not an apposition relation. This requires semantic knowledge.

We harvest some of this knowledge from the ETCBC subphrase,
which has a relation called `par`. Though the subphrases
are not very reliable, they are very useful for this particular
use-case. We run a query to find all such transitions between
substantives, and store a "separator" on the first word.
i.e. in the example above, on "Shem".
"""

import json

# this query identifies 2 words, one at the beginning of a
# subphrase, the other at the end, which are both either 
# substantives or proper nouns. The second subphrase has
# a relation of `par` meaning there is a parallel relation
# intervening between the two substantives;
# for details on the query pattern, see Text-Fabric documentation
para_query = '''

sp1:subphrase
    w1:word pdp=nmpr|subs
sp2:subphrase rela=par
    w2:word pdp#adjv|advb ls#card

sp1 := w1
sp1 <: sp2
sp1 <mother- sp2
sp2 =: w2
w1 .lex#lex. w2

'''

def get_para_seps(paths, tf_api):
    """Run subphrase query with Text-Fabric and export results."""

    # run the query and pick out the first word
    results = tf_api.S.search(para_query)
    paraseps = []
    for sp1, w1, sp2, w2 in sorted(results):
        paraseps.append(w1)

    with open(paths['paraseps'], 'w') as outfile:
        json.dump(paraseps, outfile)
