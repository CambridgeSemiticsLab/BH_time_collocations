patterns = [

# PREP + NOUN
"""
ph:phrase rela=NA
    w1:word pdp=prep
    <: w2:word pdp=subs ls#card

w1 =: ph
w2 := ph
""",

# PREP + ART + NOUN
"""
ph:phrase rela=NA
    w1:word pdp=prep
    <: word lex=H
    <: w2:word pdp=subs ls#card

w1 =: ph
w2 := ph
"""
    
# PREP + ADJV
"""
ph:phrase
    w1:word pdp=prep
    <: w2:word pdp=advb

w1 =: ph
w2 := ph
""",
    
# NOUN + C + NOUN
"""
ph:phrase rela=NA
    w1:word pdp=subs st=c ls#card
    <: w2:word pdp=subs ls#card st=a

w1 =: ph
w2 := ph
""",   
    
# PREP+ NOUN + C + NOUN
"""
ph:phrase rela=NA
    w1:word pdp=prep
    <: word pdp=subs st=c ls#card
    <: w2:word pdp=subs ls#card st=a

w1 =: ph
w2 := ph
""",   
    
# PREP + NOUN + C + NOUN + C + NOUN
"""
ph:phrase rela=NA
    w1:word pdp=prep
    <: word pdp=subs st=c ls#card
    <: w2:word pdp=subs ls#card st=a

w1 =: ph
w2 := ph
""",
    
# PREP + NOUN + C + NOUN + C + NOUN
"""
ph:phrase rela=NA
    w1:word pdp=subs st=c ls#card
    <: word pdp=subs st=c ls#card
    <: word pdp=subs st=c ls#card
    <: w2:word pdp=subs ls#card st=a

w1 =: ph
w2 := ph
""",
    
# CARD + CARD + SUBS
"""
ph:phrase rela=NA
    w1:word pdp=subs ls=card
    <: word pdp=subs ls=card
    <: w2:word pdp=subs ls#card

w1 =: ph
w2 := ph
""",
 
# CARD + W + CARD + SUBS
"""
ph:phrase rela=NA
    w1:word pdp=subs ls=card
    <: word lex=W
    <: word pdp=subs ls=card
    <: w2:word pdp=subs ls#card

w1 =: ph
w2 := ph
""",
    
# ~Cardinal quantifier phrases
"""
phrase rela=NA
/where/
    word
/have/
    pdp=subs prs=absent
/-/
    word ls=card
    word pdp=subs ls#card

""",
]