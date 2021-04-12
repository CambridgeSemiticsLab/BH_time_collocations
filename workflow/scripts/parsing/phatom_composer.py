"""
Nearly all phrase atoms have been parsed with the phatom_parser. 
These objects have relations between each other within embedding phrases
that can be accessed in BHSA through TF. The relations form a complete
(functional) phrase; whereas the atoms comprise embedded subphrases.
We will take advantage of this in order to complete the phrase parsing.
"""

import json
import collections
import copy
from tools.load_parse import ParseLoader
import tools.nav_tree as nt

class PhraseAtomComposer:

    def __init__(self, ph2parse, tf_api, protect_original=False):
        self.F, self.E, self.L = tf_api.F, tf_api.E, tf_api.L
        if protect_original:
            self.ph2parse = copy.deepcopy(ph2parse) # write only to local copy
        else:
            self.ph2parse = ph2parse
        self.mom2kids = self.build_edges()
   
    def build_edges(self):
        """Build up an edges dictionary.
        
        The ETCBC phrase atom edges are a bit complex.
        Phrase atoms can have relation edges to other 
        phrase atoms, or to individual words contained in
        another phrase atom. 

        We add edge mappings through a dict. The dict is a mapping 
        from a phrase node to either another phrase node 
        or to a tuple of slots. The mapping represents child to mother
        mapping at first. But this dict is reversed before it is 
        returned to yield a mother to children mapping.
       
        One case which we adjust here is the Link (conj) rela, 
        which points at the item that is being parallelized 
        rather that the phrase which follows the conjunction. 
        We alter that to point at the subsequent atom istead.
    
        One known issue with the atoms concerns the apposition
        relation where [NP -Appo> PP], in which the appositional 
        rela would be better pointed at a nominal element in the
        phrase. Unfortunately, the ambiguity involved with attempting
        to find the right nominal element makes solving this 
        problem beyond the scope of this project.
       """
        
        child2mom = {}
        relamap = {
            'Appo': 'APPO',
            'Spec': 'SPEC',
            'Link': 'CONJ',
            'Sfxs': 'ADJV',
            'Para': 'PARA',
            'NA': None,
        }
        
        # build the maps
        for ph in self.F.otype.s('phrase_atom'):
            
            # get data on this ph and its mother
            rela = self.F.rela.v(ph)
            rela = relamap[rela]
            mom = self.E.mother.f(ph)
            momotype = set(self.F.otype.v(m) for m in mom)
                        
            # modify Link
            if rela == 'CONJ':
                # reassign these edges to point at 
                # the parallel element instead
                child2mom[ph] = (ph+1, rela)
               
            # deal with normal phrases
            elif 'phrase_atom' in momotype:
                child2mom[ph] = (mom[0], rela)
                
            # word mothers as slots
            elif 'word' in momotype:
                child2mom[ph] = (mom, rela)
                
        # reverse the dict
        mom2kids = collections.defaultdict(list)
        for child, edge in child2mom.items():
            mom, rela = edge
            mom2kids[mom].append((child, rela))
        
        return mom2kids
        
    def get_parse(self, ph_atom):
        """Retrieve phrase atom parsing."""
        try:
            return self.ph2parse[ph_atom]
        except KeyError:
            words = self.L.d(ph_atom, 'word')
            if len(words) == 1:
                return words
            else:
                raise Exception(f'No parsing found for {ph_atom}!')
                
    def sort_children(self, node, edges):
        """Sort related nodes for walking."""
        before_mom = []
        after_mom = []
        for edge in edges:
            child, rela = edge
            if child > node:
                after_mom.append(edge)
            else:
                before_mom.append(edge)
        after_mom.sort()
        before_mom.sort(reverse=True)
        return after_mom + before_mom
        
    def compose_phrase(self, node):
        """Recursively compose phrase elements.
        
        A phrase atom in ETCBC has a set of edge relations 
        (= mother); a complete phrase is a network of these 
        relations. In this function we walk those relations 
        to retrieve the parsings that have already been produced
        by the phatom_parser which correspond with an ETCBC phrase
        atom. The parsings exist in the form of a list which contains
        3 or 1 items:

        if 3:
            [source, target, relation]
            where source and target are either BHSA slot nodes or embedded
            parse lists (=embedded custom subphrases); source has an edge 
            relation ('relation') to the target 
        elif 1:
            [integer]
            where integer is a BHSA slot node; this is a standalone phrase
            with no relations

        The network of phrase atom parsings are united by 
        modifiying the parse lists which represent them. 
        Two main kinds of modifications are made:

        1. INTERNAL modifications to parse list: thes are cases where
            another phrase atom has an edge relation to single word /
            slot in the parse; done by iteratively reassigning through 
            an index assignment operation and recursive call on a related atom
            e.g. list[1] = [recursive_call(related_atom), list[1], 'rela']
            meanwhile, the new element is built up by a
            recursive call of this function, which effectively
            replaces "new_element" with its own composed phrase, 
            and therefore drawing in any edge connections to it.

        2. EXTERNAL modifications: this is where a parse list is 
            wrapped with new relations iteratively with self-reference.
            e.g. parse = [recursive_call(related_atom), parse, 'rela']

        The result of this edge walk and composition is a parse list
        that corresponds with a complete, functional phrase according
        to the ETCBC.
        """
        
        parse = self.get_parse(node)

        # modify INTERNAL phrase constituents with slot-based edge mappings
        if len(parse) == 3:
            for ph in nt.traverse_tree(parse):                
                for i, sp in enumerate(ph[:-1]):
                    slots = tuple(sorted(nt.get_slots(sp)))
                    for kid, rela in self.mom2kids.get(slots, []):
                        ph[i] = [
                            self.compose_phrase(kid),
                            ph[i], # build up recursively
                            rela
                        ]
        # modify INTERNAL phrase constituents for single-word phrases
        elif len(parse) == 1:
            for kid, rela in self.mom2kids.get((parse[0],), []):
                parse[0] = [
                    self.compose_phrase(kid),
                    parse[0],
                    rela
                ]
        
        # single-slot parse adjustment and sanity checks
        if len(parse) == 1:
            parse = parse[0]
        elif len(parse) != 3:
            raise Exception(f'Invalid parse length of {len(parse)}: {parse}')

        # compose EXTERNAL phrase constituents by iteratively wrapping a new list
        for kid, rela in self.sort_children(node, self.mom2kids.get(node, [])):
            parse = [
                self.compose_phrase(kid),
                parse,
                rela
            ]
        
       # finish
        return parse

def get_complete_phrases(ph2parse, tf_api):
    """Retrieve phrases completely covered by the parsings.
    
    The phrase atom parser runs on phrase_atoms, which are
    component parts of a complete phrase. In some cases the
    parser was unable to parse a phrase_atom, meaning that 
    some phrases are left without a complete parsing. This
    function only selects those phrases with complete parses.
    """
    F, L = tf_api.F, tf_api.L
    parsed_atoms = set(ph2parse)

    # add unparsed conjunctions
    # these are compensated for in the Composer object
    for atom in F.otype.s('phrase_atom'):
        if F.rela.v(atom) == 'Link':
            parsed_atoms.add(atom)

    # select only those phrases completely covered by the parser
    whole_phrases = []
    for phrase in F.otype.s('phrase'):
        ph_atoms = set(L.d(phrase, 'phrase_atom'))
        if parsed_atoms.issuperset(ph_atoms):
            whole_phrases.append(phrase)
           
    return whole_phrases

def compose_phrases(paths, tf_api):
    """Compose all eligible phrases."""
    
    F, L = tf_api.F, tf_api.L
    ph2parse = ParseLoader(paths['parsed_atoms']).load()
    
    whole_phrases = get_complete_phrases(ph2parse, tf_api)
    composer = PhraseAtomComposer(ph2parse, tf_api)
    full_parses = {}
    
    # iterate through all whole phrases and call composer on the 
    # first phrase atom of each one; the network connections between
    # all of the atoms should cause them all to be grabbed
    for phrase in whole_phrases:
        first_atom = L.d(phrase, 'phrase_atom')[0]
        comp_parse = composer.compose_phrase(first_atom)

        # re-wrap single-word parses
        if type(comp_parse) == int:
            comp_parse = [comp_parse]

        # sanity check: 
        # compare slots in parse with slots in atom 
        # to make sure all slots are accounted for
        ph_slots = set(L.d(phrase,'word'))
        if len(comp_parse) == 3:
            comp_slots = set(nt.get_slots(comp_parse))
        elif len(comp_parse) == 1:
            comp_slots = set(comp_parse)
        if ph_slots != comp_slots:
            raise Exception(f'Missing slots for {phrase}; orig: {ph_slots}; comp: {comp_slots}')

        # save the parse
        full_parses[phrase] = comp_parse
            
    # export
    with open(paths['parsed_phrases'], 'w') as outfile:
        json.dump(full_parses, outfile, indent=2)

#    whole_atoms= set(at for ph in whole_phrases for at in L.d(ph,'phrase_atom') if at in parsed_atoms)
#    part_atoms = parsed_atoms - whole_atoms
