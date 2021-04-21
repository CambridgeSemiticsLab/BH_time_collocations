"""
Provide functions for tracking and saving
annotations that have been accepted already.
"""

import json
import dictdiffer
from IPython.display import HTML, display
from pathlib import Path

def load_style():
    """Load styles for BHSA divs."""
    display(HTML(Path('bhsa.css').read_text()))

class Tracker:

    def __init__(self, origpath, trackpath, silent=False):
        self.trackpath = trackpath
        self.origpath = origpath
        self.silent = silent
        with open(trackpath, 'r') as infile:
            self.tracked = json.load(infile)
        with open(origpath, 'r') as infile:
            self.orig = json.load(infile)

    def save(self, clause, parse):
        """Save and keep track of a given parse."""
        if type(clause) == int:
            clause = str(clause)

        orig_parse = self.orig[clause]
        corrections = list(dictdiffer.diff(orig_parse, parse))
        self.tracked[clause] = corrections

        with open(self.trackpath, 'w') as outfile:
            json.dump(self.tracked, outfile, indent=2)
        if not self.silent:
            spans = [
                f'clause {clause} saved (total {len(self.tracked)})'
            ]
            if corrections:
                spans.extend(str(c) for c in corrections)
            else:
                spans.append('no corrections')
            spans = '<br>'.join(spans)
            display(HTML(f'<div style="background:#9BD788;display:inline-block">{spans}</div>'))
