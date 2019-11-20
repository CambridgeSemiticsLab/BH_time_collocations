import sys

class Debugger(object):
    """Display debugging messages if toggled"""
    def __init__(self, boolean):
        self.report = boolean
        self.indent = 0
    def say(self, msg, indent=0, end='\n'):
        self.indent = indent or self.indent
        if self.report:
            indent = self.indent * '\t'
            fmtmsg = f'{indent}{msg}{end}'
            sys.stderr.write(fmtmsg)