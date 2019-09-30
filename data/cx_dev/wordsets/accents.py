import re
import collections
from positions import Positions

def book_class(node, tf):
    '''
    Returns the accent class of a node's
    book, i.e. the 21 or the 3
    '''
    book = tf.api.T.sectionFromNode(node)[0]
    if book not in ('Psalms', 'Job', 'Proverbs'):
        return '21'
    else:
        return '3'
    
def masoretic_word(word, tf):
    '''
    Retrieves complete phonological unit
    (thanks to Johan Lundberg for terminology here).
    Returns sorted list of BHSA word nodes.
    
    If a word is followed by maqqeph (ETCBC "&")
    or zero-space it is part of a larger phonological unit.
    But if a word has maqqeph and its own accent, it
    is treated separately from the subsequent word.
    '''
    F, T = tf.api.F, tf.api.T
    
    def maketext(word):
        # make text for regexing
        return str(F.trailer.v(word))

    mwords = {word} # collect them here
    thisword = word-1
    text = maketext(thisword)

    # back up `this_word` to beginning
    while ('&' in text) or (F.trailer.v(thisword) == ''):
        mwords.add(thisword)
        thisword = thisword-1
        text = maketext(thisword)

    # restart at middle
    thisword = word
    text = maketext(thisword)

    # move from middle to end
    while ('&' in text) or (F.trailer.v(thisword) == ''):
        mwords.add(thisword+1)
        thisword = thisword+1
        text = maketext(thisword)

    return tuple(sorted(mwords))

class Accents:
    
    '''
    Classifies words as disjunct or conjunct
    depending on the cantillation signs.
    Makes these sets available under
    a series of dict mappings, allowing the
    results to be scrutinized.
    '''
    
    def __init__(self, tf):
        # set up Text-Fabric classes
        self.tf = tf
        self.F, self.T, self.L = tf.api.F, tf.api.T, tf.api.L
        
        # disjunctive accent patterns (ETCBC transcription)
        disA = {
            '21': {
                'paseq': '.*05',
                'atnach': '.*92',
                'silluq': '.*75',
                #'tiphchah' see special function below
                'zaqeph qaton': '.*80',
                'zaqeph gadol': '.*85',
                'segolta': '.*01',
                'shalshelet': '.*65',
                'rebia': '.*81',
                'zarqa': '.*02',
                'pashta': '.*03',
                'yetiv': '.*10', 
                'tebir': '.*91',
                'geresh': '.*(61|11)',
                'gershayim': '.*62',
                'pazer qaton': '.*83',
                'qarney parah': '.*84',
                'telisha gedola': '.*(14|44)',
            },
            '3': {
                'paseq': '.*05',
                'atnach': '.*92',
                'silluq': '.*75',
                'rebia': '.*81',
                'oleh weyored': '.*60.*71',
                'rebia mugrash': '.*11.*81',
                'shalshelet gedolah': '.*65.*05',
                'tsinor': '.*82',
                'dechi': '.*13',
                'pazer':  '.*83',
                'mehuppak legarmeh': '.*70.*05',
                'azla legarmeh': '.*(63|33).*05'
            }
        }
        
        # conjunct accent patterns
        conA = {
            '21': {
                'munach': '.*74',
                'mehuppak': '.*70',
                'mereka': '.*71',
                'merekah kefula': '.*72',
                'darga': '.*94',
                'azla/qadma': '.*(63|33)',
                'telisha qetannah': '.*04',
                'yerah': '.*93',
                'mayela': '.*73\S+(75|92)', # assumes _ replaced with \s
            },
            '3': {
                'munach': '.*74',
                'mereka': '.*71',
                'illuy': '.*64',
                'tarcha (tiphcha)': '.*73',
                'yerah': '.*93',
                'mehuppak': '.*70',
                'azla/qadma': '.*(63|33)',
            }
        }
        
        # compile dis regex pattern
        self.disRE = {bclass: {name:re.compile(patt) for name, patt in names.items()} 
                          for bclass, names in disA.items()}
        self.tiphchah_RE = re.compile('.*73') # requires special check
        
        # compile con regex pattern
        self.conRE = {bclass: {name:re.compile(patt) for name, patt in names.items()} 
                          for bclass, names in conA.items()}
        
        # generate masoretic word sets
        self.mwords = {}
        for w in self.F.otype.s('word'):
            self.mwords[w] = masoretic_word(w, self.tf)
        
        # assemble word sets
        self.accenttype = {}
        self.atype2set = collections.defaultdict(list)
        self.atype2name2set = collections.defaultdict(lambda:collections.defaultdict(list))
        
        # loop and assign labels
        for w in set(self.F.otype.s('word')):
            dismatches = self.disjunct(w)
            if not dismatches: # only do if necessary
                conmatches = self.conjunct(w)                
            if dismatches:
                self.accenttype[w] = 'disjunct'
                self.atype2set['disjunct'].append((w,))
                self.atype2name2set['disjunct'][dismatches].append((w,))
            elif conmatches:
                self.accenttype[w] = 'conjunct'
                self.atype2set['conjunct'].append((w,))
                self.atype2name2set['conjunct'][conmatches].append((w,))
            else:
                self.accenttype[w] = 'unknown'
                self.atype2set['unknown'].append((w,))
                
    def clean(self, text):
        '''
        Replaces certain transcriptions.
        '''
        return text.replace('_', ' ')
        
    def disjunct(self, word):
        '''
        Evaluates simple cases of disjunction
        with a regex match.
        '''
        # get and test phonological unit
        mword = self.T.text(self.mwords[word], fmt='text-trans-full')
        mword = self.clean(mword)
        bclass = book_class(word, self.tf)
        
        # identify and return matches 
        matches = []
        for name, patt in self.disRE[bclass].items():
            if patt.match(mword):
                matches.append(name)
        if self.tiphchah(mword):
            matches.append('tiphchah')
            
        return tuple(matches)
    
    def tiphchah(self, mword):
        '''
        There are rare cases where a tiphchah
        is re-evaluated as a mayelah, i.e. 
        when tiphchah is on a word with atnach or silluq.
        This function simply checks for that case and 
        excludes it in order to validate tiphchah.
        '''
        if (self.tiphchah_RE.match(mword) 
                and not self.conRE['21']['mayela'].match(mword)):
            return True
    
    def conjunct(self, word):
        '''
        Returns a list of conjunctive accent matches.
        !!CAUTION!! Should only be used after a negative
        test for disjunctive accents. Many conjunctive
        accents belong to larger patterns that must first
        be checked. In this class, this method is used
        in an `elif` only AFTER checking disjunctives.
        '''
        # get and test phonological unit
        mword = self.T.text(self.mwords[word], fmt='text-trans-full')
        mword = self.clean(mword)
        bclass = book_class(word, self.tf)
        
        # identify and return matches 
        matches = []
        for name, patt in self.conRE[bclass].items():
            if patt.match(mword):
                matches.append(name)
        return tuple(matches)
        

dis21 = {
    'atnach': {'regex':'.*\u0591', 'etcbc':'92'},
    'tiphchah': {'regex':'.*\u0596', 'etcbc':'73'},
    'zaqeph qaton': {'regex':'.*\u0594', 'etcbc':'80'},
    'zaqeph gadol': {'regex':'.*\u0595', 'etcbc':'85'},
    'segolta': {'regex':'.*\u0592', 'etcbc':'01'},
    'shalshelet': {'regex':'.*\u0593', 'etcbc':'65'},
    'rebia': {'regex':'.*\u0597', 'etcbc':'81'},
    'zarqa': {'regex':'.*\u05AE', 'etcbc':'02'},
    'pashta': {'regex':'.*\u0599', 'etcbc':'03'},
    'yetiv': {'regex':'.*\u059A', 'etcbc':'10'}, 
    'tebir': {'regex':'.*\u059B', 'etcbc':'91'},
    'geresh': {'regex':'.*\u059C', 'etcbc':'61'},
    'gershayim': {'regex':'.*\u059E', 'etcbc':'62'},
    'legarmeh': {'regex':'.*\u05A3.*\u05C0', 'etcbc':'[74...05]'},
    'pazer qaton': {'regex':'.*\u05A1', 'etcbc':'83'},
    'pazer gadol': {'regex':'.*\u059F', 'etcbc':'84'},
    'telisha gedola': {'regex':'.*\u05A0', 'etcbc':'14|44'},
}
dis3 = {
    'atnach': {'regex':'.*\u0591', 'etcbc':'92'},
    'rebia': {'regex':'.*\u0597', 'etcbc':'81'},
    'oleh weyored': {'regex':'.*\u05AB.*\u05A5', 'etcbc':'.*60.*71'},
    'rebia mugrash': {'regex':'.*\u059D.*\u0597', 'etcbc':'.*11.*81'},
    'shalshelet gedolah': {'regex':'.*\u0593.*\u05C0', 'etcbc':'[65...05]'},
    'tsinor': {'regex':'.*\u0598', 'etcbc':'82'},
    'dechi': {'regex':'.*\u05AD', 'etcbc':'13'},
    'pazer': {'regex':'.*\u05A1', 'etcbc':'83'},
    'mehuppak legarmeh': {'regex':'.*\u05A4.*\u05C0', 'etcbc':'[70...05]'},
    'azla legarmeh': {'regex':'.*\u05A8.*\u05C0', 'etcbc':'[63|33...05]'}
}

# conjunct UTF8 patterns
conA = {
    '21': {
        'munach': '.*\u05A3',
        'mehuppak': '.*\u05A4',
        'mereka': '.*\u05A5',
        'darga': '.*\u05A7',
        'azla/qadma': '.*\u05A8',
        'telisha qetannah': '.*\u05A9',
        'yerah': '.*\u05AA',
    },
    '3': {
        'munach': '.*\u05A3',
        'mereka': '.*\u05A5',
        'illuy': '.*\u05AC',
        'tarcha': '.*\u0596',
        'yerah': '.*\u05AA',
        'mehuppak': '.*\u05A4',
        'azla/qadma': '.*\u05A8',
    }
}