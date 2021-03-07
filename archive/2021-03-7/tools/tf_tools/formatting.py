'''
Functions that perform
various helper tasks with Text-Fabric
'''

import collections

# -- SBL-style References Formatting --

book2sbl = {
'Genesis': 'Gen', 'Exodus': 'Exod', 'Leviticus': 'Lev', 'Numbers': 'Num',
'Deuteronomy': 'Deut', 'Joshua': 'Josh', 'Judges': 'Judg', '1_Samuel': '1 Sam', '2_Samuel': '2 Sam',
'1_Kings': '1 Kgs', '2_Kings': '2 Kgs', 'Isaiah': 'Isa', 'Jeremiah': 'Jer', 'Ezekiel': 'Ezek',
'Hosea': 'Hos', 'Joel': 'Joel', 'Amos': 'Amos', 'Obadiah': 'Obad', 'Jonah': 'Jonah', 'Micah': 'Mic',
'Nahum': 'Nah', 'Habakkuk': 'Hab', 'Zephaniah': 'Zeph', 'Haggai': 'Hag', 'Zechariah': 'Zech',
'Malachi': 'Mal', 'Psalms': 'Ps', 'Job': 'Job', 'Proverbs': 'Prov', 'Ruth': 'Ruth',
'Song_of_songs': 'Song', 'Ecclesiastes': 'Eccl', 'Lamentations': 'Lam', 'Esther': 'Esth',
'Daniel': 'Dan', 'Ezra': 'Ezra', 'Nehemiah': 'Neh', '1_Chronicles': '1 Chr', '2_Chronicles': '2 Chr'}

def formatPassages(resultslist, tf):
    '''
    Formats biblical passages with SBL style
    for a list of results.
    '''
    
    T = tf.T # get TF Text Class
    
    # compile results in nested dict
    # --enables formatting to be on a 
    # book and chapter by chapter basis
    book2ch2vs = collections.defaultdict(lambda: collections.defaultdict(set))
    for result in sorted(resultslist):
        book, chapter, verse = T.sectionFromNode(result[0])
        book = book2sbl[book]                
        book2ch2vs[book][chapter].add(str(verse))
            
    # assemble in to readable passages list
    passages = []
    for book, chapters in book2ch2vs.items():
        ch_verses = []
        for chapter, verses in chapters.items():
            verses = ', '.join(f'{chapter}:{verse}' for verse in sorted(verses))
            ch_verses.append(verses)
        passage = f'{book} {", ".join(ch_verses)}'
        passages.append(passage)
            
    return '; '.join(passages)

# -- Pretty Function Formatting -- 
# converts ETCBC codes to long-style function strings

funct2function = '''
Adju	Adjunct	
Cmpl	Complement	
Conj	Conjunction	
EPPr	Enclitic personal pronoun	
ExsS	Existence with subject suffix	
Exst	Existence	
Frnt	Fronted element	
Intj	Interjection	
IntS	Interjection with subject suffix	
Loca	Location	
Modi	Modifier	
ModS	Modifier with subject suffix	
NCop	Negative copula	
NCoS	Negative copula with subject suffix	
Nega	Negation	
Objc	Object	
PrAd	Predicative adjunct	
PrcS	Predicate complement with subject suffix	
PreC	Predicate complement	
Pred	Predicate	
PreO	Predicate with object suffix	
PreS	Predicate with subject suffix	
PtcO	Participle with object suffix	
Ques	Question	
Rela	Relative	
Subj	Subject	
Supp	Supplementary constituent	
Time	Time	
Unkn	Unknown	
Voct	Vocative
'''.split('\n')

funct2function = dict((funct[0], funct[1].strip('\t')) for funct in [func.split('\t', 1) for func in funct2function] if funct[0])

# simplified function labels
simplified_functions = {
    'PreO': 'Pred', 
    'PreS': 'Pred', 
    'PtcO': 'Pred',
    'IntS': 'Intj', 
    'NCoS': 'NCop',
    'ModS': 'Modi',
    'ExsS': 'Exst'
}

