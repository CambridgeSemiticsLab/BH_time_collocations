import re
import collections

class AccentTagger:
    
    """Classifies word accent as disjunct or conjunct.

    This is done by looking up the word's transcription
    and then looking for regex matches for particular 
    codes in the transcribed text.

    The source_data comes from the ETCBC's BHSA and is processed
    with Text-Fabric methods. The transcription is the
    ETCBC's, which includes codes for the various accents.
    The codes are matched with a series of regular expressions.

    Run tag on a BHSA word node number to get a string:
        conjunct
        disjunct
        unknown (to cover edge cases)
    """
    
    def __init__(self, tf_api):

        # set up Text-Fabric classes
        self.tf_api = tf_api
        self.F, self.T, self.L = tf_api.F, tf_api.T, tf_api.L
        
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
       
    def masoretic_word(self, word):
        """
        Retrieves complete phonological unit
        (thanks to Johan Lundberg for terminology here).
        Returns sorted list of BHSA word nodes.
        
        If a word is followed by maqqeph (ETCBC "&")
        or zero-space it is part of a larger phonological unit.
        But if a word has maqqeph and its own accent, it
        is treated separately from the subsequent word.
        """
        F, T = self.F, self.T
        
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
               
    def clean(self, text):
        """
        Replaces certain transcriptions.
        """
        return text.replace('_', ' ')

    def transcribe_word(self, word, **kwargs):
        """Transcribe a word in its phonological context."""

        # get the transcription of the 
        # word's phonological unit
        if kwargs.get('phono', True):
            word = self.masoretic_word(word)
        else:
            word = [word]
            
        mword_trans = self.T.text(
            word,
            fmt='text-trans-full'
        )
        mword_trans = self.clean(mword_trans)
        return mword_trans

    def book_class(self, node):
        """
        Returns the accent class of a node's
        book, i.e. the 21 or the 3
        """
        book = self.T.sectionFromNode(node)[0]
        if book not in ('Psalms', 'Job', 'Proverbs'):
            return '21'
        else:
            return '3'
 
    def tiphchah(self, mword):
        """
        There are rare cases where a tiphchah
        is re-evaluated as a mayelah, i.e. 
        when tiphchah is on a word with atnach or silluq.
        This function simply checks for that case and 
        excludes it in order to validate tiphchah.
        """
        if (self.tiphchah_RE.match(mword) 
                and not self.conRE['21']['mayela'].match(mword)):
            return True
 
    def disjunct(self, word, **kwargs):
        """
        Evaluates simple cases of disjunction
        with a regex match.
        """
        mword_trans = self.transcribe_word(word, **kwargs)
        bookclass = self.book_class(word)

        # identify matches
        matches = []
        dispatterns = self.disRE
        for name, patt in dispatterns[bookclass].items():
            if patt.match(mword_trans):
                matches.append(name)
        if self.tiphchah(mword_trans):
            matches.append('tiphchah')
            
        return matches
   
    def conjunct(self, word, **kwargs):
        """
        Returns a list of conjunctive accent matches.
        !!CAUTION!! Should only be used after a negative
        test for disjunctive accents. Many conjunctive
        accents belong to larger patterns that must first
        be checked. In this class, this method is used
        in an `elif` only AFTER checking disjunctives.
        """
        mword_trans = self.transcribe_word(word, **kwargs)
        bookclass = self.book_class(word)
        
        # identify matches 
        matches = []
        conpatterns = self.conRE
        for name, patt in conpatterns[bookclass].items():
            if patt.match(mword_trans):
                matches.append(name)

        return matches

    def tag(self, w, **kwargs):
        """Tag a word's accent in its phonological context.

        Returns a string tag.
        """
       
        # look for conjunctive/disjunctive accents
        # we only check for conjunctive accents
        # if no disjunctive ones are found
        dismatches = self.disjunct(w, **kwargs)
        if not dismatches:
            conmatches = self.conjunct(w, **kwargs)      

        # return the tag
        if dismatches:
            return 'disjunct'
        elif conmatches:
            return 'conjunct'
        else:
            return 'unknown'
