import re
import collections
from itertools import groupby, count

def as_range(iterable):
    l = list(iterable)
    if len(l) > 1:
        return '{0}-{1}'.format(l[0], l[-1])
    else:
        return '{0}'.format(l[0])

def ranges(numbers, joiner=', '):
    groups = groupby(numbers, lambda n, c=count(): n-next(c))
    return joiner.join(
        as_range(g) for _, g in groups
    )

def get_verserefs(verses):
    """Compile verse-ref string."""
    
    # ensure book insertion order with Ordered Dict
    bk2ch2vs = collections.OrderedDict()
    
    # make first pass to cluster verses
    for verse in verses:
        bk, ch, vs = re.match('^(.*) (\d+):(\d+)', verse).groups()
        if bk not in bk2ch2vs:
            bk2ch2vs[bk] = collections.OrderedDict()
        if ch not in bk2ch2vs[bk]:
            bk2ch2vs[bk][ch] = set()
        bk2ch2vs[bk][ch].add(vs)
        
    # compile string
    bkstrs = []
    for bk,chs in bk2ch2vs.items():
        chstrs = []
        for ch, vss in chs.items():
            # compress consecutive verses to range
            vss = ranges(sorted(int(i) for i in vss))
            chstrs.append(ch+':'+vss)
        bkstrs.append(bk+' '+'; '.join(chstrs))
    
    refstring = '; '.join(bkstrs)
    
    return refstring