"""
Information display with timestamps
"""

from datetime import datetime

class Info:
    """Provide timed messages"""
    def __init__(self):
        self.start = datetime.now()
        self.indent = 0
    def timestamp(self):
        # give elapsed time 
        return f'{datetime.now()-self.start}'
    def msg(self, m, indent=0):
        indent = indent or self.indent
        indent = ' | \t' * indent
        print(f'{indent}{timestamp()} {m}')
