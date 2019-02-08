'''
This module contains parameters to be used
with contextcount.ContextCounter. The class DeliverParameters
simply allows the parameters to be configured with a 
supplied list of target terms.

The parameters are packaged into dictionaries and
delivered under the variable of 'parameters' in the form 
of an iterable of all the dictionaries. 
Each dictionary contains all the parameters needed to query and process a particular syntactic context.

The parameters consist of a search string (template),
tokenizer functions for indexed results from running the
search templates, and optional arguments for formatting
or analyzing the templates. 

The following keys are expected for each parameter:

template - a TF search template (string)
target - the index of the target word
bases - tuple of basis indexes
target_tokenizer - a function to construct target tokens, requires index of target
basis_tokenizer - a function to construct basis tokens, requires index of basis
*kwargs - optional dict of keyword arguments for formatting the template (if any)
*sets - optional dict containing sets arguments for TF.search
*collapse_instances - optional boolean on whether to collapse multiple instances of a 
                      basis element at the result level; default is False

'''

from textwrap import dedent # removes identations from triple-quote strings to align them with function
import re

def deliver_params(target_nouns, tf=None):
    '''
    A function that simply
    umbrellas all of the parameters
    and allows a target noun list to 
    be fed to the params without import
    statements. Returns a list of dictionaries
    of parameters. TF is an instance of Text-Fabric.
    '''
    F, L = tf.F, tf.L # text-fabric feature class
    parameters = list()
    sets = {'target': target_nouns}

    # <><> universal tokenizers (i.e. used throughout) <><>

    def disambigUTF8(wordnode):
        '''
        Adds modified ETCBC disambiguators 
        to UTF8 lexical forms.
        '''
        letters = [letter for letter in F.lex.v(wordnode)]
        gloss_number = letters.count('=') + 1 if {'[', '/'}&set(letters) else '' # count = (etcbc disambiguator)
        wordtype = 'n' if '/' in letters else 'v' if '[' in letters else ''
        utf8 = F.lex_utf8.v(wordnode)
        disambigs = [dis for dis in [utf8, wordtype+str(gloss_number)] if dis] # remove null strings (e.g. in cases of preps)
        return '.'.join(disambigs)
    
    def token_lex(target):
        '''
        Builds simple lexeme token 
        with disambiguation if necessary.
        '''
        lexeme = disambigUTF8(target)
        return lexeme

    def token_verb(verb_node):
        '''
        Constructs a stem+lemma verb token.
        '''
        verb_stem = F.vs.v(verb_node)
        verb_lex = disambigUTF8(verb_node)
        return f'{verb_lex}.{verb_stem}'


    # <><> clause constituent relation searches <><>


    # target in a non-prepositional function relating to a predicate lexeme
    clause_function = '''

    clause kind=VC
    /without/
        phrase function=PreC
    /-/
        phrase function=Pred|PreO|PreS
            <head- word pdp=verb lex#HJH[
        phrase
        /with/
        typ#PP
        /or/
        typ=PP function=Objc
        /-/
            <nhead- target

    % target=4
    % bases=(2, 3)
    '''

    def token_cc_rela(bases, target):
        # T.function→ st.verb.lex
        target_funct = F.function.v(bases[1])
        verb = token_verb(bases[0])
        return f'T.{target_funct}→ {verb}'

    parameters.append({'template': dedent(clause_function), 
                       'target': 4, 
                       'bases': (2,3), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_cc_rela,
                       'sets': sets,
                       'name': 'T.function→ st.verb.lex'
                      })


    # target in a prepositional function relating to a predicate lexeme
    clause_function_PP = '''

    clause kind=VC
    /without/
        phrase function=PreC
    /-/
        phrase function=Pred|PreO|PreS
            <head- word pdp=verb lex#HJH[

        phrase typ=PP function#Objc
            <head- word
            <obj_prep- target

    % target=5
    % bases=(2,3,4)
    '''

    def token_cc_rela_PP(bases, target):
        # T.prep.funct→ st.verb.lex
        verb, prep = token_verb(bases[0]), token_lex(bases[2])
        target_funct = F.function.v(bases[1])    
        return f'T.{prep}.{target_funct}→ {verb}'

    parameters.append({'template': dedent(clause_function_PP), 
                       'target': 5, 
                       'bases': (2,3,4), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_cc_rela_PP,
                       'sets': sets,
                       'name': 'T.prep.funct→ st.verb.lex'
                      })
    
    
    # target is subject to a complement lexeme in a predicate-complement relation 
    ccr_subj_preC = '''

    clause
        phrase function=Subj
            <head- target
        phrase function=PreC typ#PP
            <head- word pdp=subs|nmpr|verb|adjv|advb sem_set#quant

    % target=2
    % bases=(4,)
    '''

    def token_subj_preC(bases, target):
        # lex.PreC→ T.Subj
        compliment = token_lex(bases[0])
        return f'{compliment}.PreC→ T.Subj'

    parameters.append({'template': dedent(ccr_subj_preC), 
                       'target': 2, 
                       'bases': (4,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_subj_preC,
                       'sets': sets,
                       'name': 'lex.PreC→ T.Subj'
                      })
    
    
    # target is subject to a prepositional complement lexeme in a predicate-complement relation
    ccr_subj_preC_PP = '''

    clause
        phrase function=Subj
            <head- target
        phrase function=PreC typ=PP
            <head- word
            <obj_prep- word pdp=subs|nmpr|verb|adjv|advb sem_set#quant

    % target=2
    % bases=(4, 5)
    '''

    def token_subj_preC_PP(bases, target):
        # lex.prep.PreC→ T.Subj
        prep, complement = token_lex(bases[0]), token_lex(bases[1])
        return f'{complement}.{prep}.PreC→ T.Subj'

    parameters.append({'template': dedent(ccr_subj_preC_PP), 
                       'target': 2, 
                       'bases': (4,5), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_subj_preC_PP,
                       'sets': sets,
                       'name': 'lex.prep.PreC→ T.Subj'
                      })

    
    # target is complement to a subject lexeme in predicate-complement relation
    ccr_preC_subj = '''

    clause
        phrase function=Subj
            <head- word pdp=subs|nmpr|verb|adjv|advb sem_set#quant
        phrase function=PreC typ#PP
            <head- target

    % target=4
    % bases=(2,)
    '''

    def token_PreC_subj(bases, target):
        # T.PreC→ Subj.lex
        subj = token_lex(bases[0])
        return f'T.PreC→ {subj}.Subj'

    parameters.append({'template': dedent(ccr_preC_subj), 
                       'target': 4, 
                       'bases': (2,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_PreC_subj,
                       'sets': sets,
                       'name': 'T.PreC→ lex.Subj'
                     })

    
    # target is a prepositional complement to a subject lexeme in predicate-complement relation
    ccr_preC_subj_PP = '''

    clause
        phrase function=Subj
            <head- word
        phrase function=PreC typ=PP
            <head- word
            <obj_prep- target

    % target=5
    % bases=(2,4)
    '''

    def token_PreC_subj_PP(bases, target):
        # T.prep.PreC→ lex.Subj
        subj, prep = token_lex(bases[0]), token_lex(bases[1])
        return f'T.{prep}.PreC→ {subj}.Subj'

    parameters.append({'template': dedent(ccr_preC_subj_PP), 
                       'target': 5, 
                       'bases': (2,4), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_PreC_subj_PP,
                       'sets': sets,
                       'name': 'T.prep.PreC→ lex.Subj'
                     })

    # <><> subphrase relation searches <><>


    #  a lexeme connected to target through a conjunction
    sp_parallel = '''

    clause
        p1:phrase
            subphrase
                w1:target
            <mother- subphrase rela=par
                w2:word sem_set#quant


    p1 <nhead- w1
    p1 <nhead- w2

    % target=3
    % bases=(5,)
    '''

    def token_sp_parallel(bases, target):
        # lex.coord→ T
        par = token_lex(bases[0])
        return f'{par}.coord→ T'

    parameters.append({'template': dedent(sp_parallel), 
                       'target': 3, 
                       'bases': (5,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_sp_parallel,
                       'sets': sets,
                       'name': 'lex.coord→ T'
                      })


    #  target connected to a lexeme through a conjunction
    sp_parallel_rela = '''

    clause
        p1:phrase
            subphrase
                w1:word sem_set#quant
            <mother- subphrase rela=par
                w2:target


    p1 <nhead- w1
    p1 <nhead- w2

    % target=5
    % bases=(3,)
    '''

    def token_sp_parallel_rela(bases, target):
        # T.coord→ lex
        par = token_lex(bases[0])
        return f'T.coord→ {par}'

    parameters.append({'template': dedent(sp_parallel_rela), 
                       'target': 5, 
                       'bases': (3,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_sp_parallel_rela,
                       'sets': sets,
                       'name': 'T.coord→ lex'
                      })


    # target is described adjectivally by a lexeme 
    sp_adjective = '''

    clause
        p1:phrase
            subphrase
                w1:target
            <mother- subphrase rela=atr
            /without/
                word ls=card|ordn
            /-/
                word pdp=adjv lex#KL/

    p1 <nhead- w1

    % target=3
    % bases=(5,)
    '''

    def token_sp_adjective(bases, target):
        # lex.atr→ T
        adjv = token_lex(bases[0])
        return f'{adjv}.atr→ T'

    parameters.append({'template': dedent(sp_adjective), 
                       'target': 3, 
                       'bases': (5,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_sp_adjective,
                       'sets': sets,
                       'name': 'lex.atr→ T'
                      })

    
    '''
    **NOTES REGARDING CONSTRUCT RELATIONS**

    !! CONSTRUCT RELATIONS HAVE BEEN COMMENTED OUT FROM THE ANALYSIS !!

    Running the model with construct relations has shown
    that they introduce a lot of noisy information. For instance,
    compare the relatedness of a construction like "sons and fathers," 
    wherein both "son" and "father" are obviously related, with the 
    less-informative "man of war" and "man of YHWH" (example from the
    previous iteration of data). This construction, in particular, 
    resulted in YHWH being rated similar to war (מלחמה). The ultimate
    issue with using the construct relation for semantic extraction is
    its ambiguity. The construct can subtly encode a whole host of 
    different kinds of relations. In "man of war," the man is portrayed
    perhaps as a subject of an implied verb "war maker," or some similar
    relation. In "man of YHWH," the "genitive" is more of a possessive 
    genitive.

    There is hope for reintegrating the construct relations into the analysis.
    For instance, after analyzing YHWH and מלחמה, it might be possible to abstract 
    semantic scores or classes that can be applied in place of the lexemes. A
    formal characteristic of the possessive use might be found in the 
    person-ness of the genitive noun. Future research might thus be able 
    to cluster various construct uses, and the classes can in turn be used to
    further refine the overall model.
    '''
    
    
    
    # a lexeme is construct with target
    sp_construct = '''

    clause
        phrase
            target st=c
            {article}
            <: basis:word pdp=subs sem_set#quant

    % target=2
    % bases=(3|4,)
    '''

    def token_sp_construct(bases, target):
        # lex.const→ T
        constr = token_lex(bases[0])
        return f'{constr}.const→ T'

#     parameters.append({'template': dedent(sp_construct), 
#                        'target': 2, 
#                        'bases': (3,), 
#                        'target_tokenizer': token_lex, 
#                        'basis_tokenizer': token_sp_construct,
#                        'sets': sets,
#                        'kwargs': {'article':''},
#                        'name': 'lex.const→ T'
#                       })

#     # with article separation
#     parameters.append({'template': dedent(sp_construct), 
#                        'target': 2, 
#                        'bases': (4,), 
#                        'target_tokenizer': token_lex, 
#                        'basis_tokenizer': token_sp_construct,
#                        'sets': sets,
#                        'kwargs': {'article':'<: word pdp=art'},
#                        'name': 'lex.const→ T (with article separation)'
#                       }) 
    
    # target is construct with a lexeme
    sp_construct_rela = '''

    clause
        phrase
            word sem_set#quant|prep pdp#art st=c
            {article}
            <: target
                    
    % target=3|4
    % bases=(2,)
    '''

    def token_sp_construct_rela(bases, target):
        # T.const→ lex
        absolute = token_lex(bases[0])
        return f'T.const→ {absolute}'
    
#     parameters.append({'template': dedent(sp_construct_rela), 
#                        'target': 3, 
#                        'bases': (2,), 
#                        'target_tokenizer': token_lex, 
#                        'basis_tokenizer': token_sp_construct_rela,
#                        'kwargs': {'article':''},
#                        'sets': sets,
#                        'name': 'T.const→ lex'
#                       })
    
#     parameters.append({'template': dedent(sp_construct_rela), 
#                        'target': 4, 
#                        'bases': (2,), 
#                        'target_tokenizer': token_lex, 
#                        'basis_tokenizer': token_sp_construct_rela,
#                        'kwargs': {'article':'<: word pdp=art'},
#                        'sets': sets,
#                        'name': 'T.const→ lex (with article separation)'
#                       })


    # <><> phrase atom relation searches <><>


    # a lexeme is connected through conjunction to target (in another phrase_atom)
    pa_parallel = '''

    clause
        p:phrase
            phrase_atom
                w1:target
            <mother- phrase_atom rela=Para
                w2:word pdp=subs sem_set#quant
    
    p <nhead- w1
    p <nhead- w2
    
    % target=3
    % bases=(5,)
    '''
    
    def token_pa_parallel(bases, target):
        # lex.coord→ T
        para = token_lex(bases[0])
        return f'{para}.coord→ T'

    parameters.append({'template': dedent(pa_parallel), 
                       'target': 3, 
                       'bases': (5,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_pa_parallel,
                       'sets': sets,
                       'name': 'lex.coord→ T (phrase atoms)'
                      })


    # target is connected through conjunction to a lexeme
    pa_parallel_rela = '''

    clause
        p:phrase
            phrase_atom
                w1:word pdp=subs sem_set#quant
            <mother- phrase_atom rela=Para
                w2:target

    p <nhead- w1
    p <nhead- w2

    % target=5
    % bases=(3,)
    '''

    def token_pa_parallel_rela(bases, target):
        # T.coord→ lex
        paralleled = token_lex(bases[0])
        return f'T.coord→ {paralleled}'

    parameters.append({'template': dedent(pa_parallel_rela), 
                       'target': 5, 
                       'bases': (3,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_pa_parallel_rela,
                       'sets': sets,
                       'name': 'T.coord→ lex (phrase atoms)'
                      })


    # a lexeme is in apposition to target word
    # NB that for appositions I opt for the older noun_heads feature,
    # since new heads does not yet include phrase_atom heads, and those
    # are a necessity for apposition relations
    pa_apposition = '''

    clause
        phrase_atom
            -noun_heads> target
        <mother- phrase_atom rela=Appo
            -noun_heads> word pdp=subs|nmpr sem_set#quant
                
    % target=2
    % bases=(4,)
    '''

    def token_pa_apposition(bases, target):
        # lex.appo→ T
        appo = token_lex(bases[0])
        return f'{appo}.appo→ T'

    parameters.append({'template': dedent(pa_apposition), 
                       'target': 2, 
                       'bases': (4,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_pa_apposition,
                       'sets': sets,
                       'name': 'lex.appo→ T'
                      })


    # target is in apposition to a lexeme
    pa_apposition_rela = '''

    clause
        phrase_atom
            -noun_heads> word sem_set#quant
        <mother- phrase_atom rela=Appo
            -noun_heads> target

    % target=4
    % bases=(2,)
    '''

    def token_pa_apposition_rela(bases, token):
        # T.appo→ lex
        appo_rela = token_lex(bases[0])
        return f'T.appo→ {appo_rela}'

    parameters.append({'template': dedent(pa_apposition_rela), 
                       'target': 4, 
                       'bases': (2,), 
                       'target_tokenizer': token_lex, 
                       'basis_tokenizer': token_pa_apposition_rela,
                       'sets': sets,
                       'name': 'T.appo→ lex'
                      })
    
    return parameters
