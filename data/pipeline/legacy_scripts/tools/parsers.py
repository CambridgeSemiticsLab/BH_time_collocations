import json
from .positions import PositionsTF

class PositionsParser:
    """Parse items from a Text-Fabric resource.

    This class provides methods for parsing basic tags
    based on contextual source_data. The contexts are accessed
    through a custom-made Positions class, which allows for
    positions (node) before or after a given node to be 
    looked up and referred to.

    Parsings are defined by creating a new class that
    inherits this one. Methods are written in the new
    class that return a boolean on whether a given item
    is a match for a given parse. The identity of that 
    parse is determined from the method's name.

    NOTE: all methods and attributes without a prefixed 
    underscore will be treated as a parsing function. To
    avoid errors, ensure that any custom methods added
    to the class which are not also a parser are prefixed
    with an underscore.
    """

    def __init__(self, tf_api):

        """Initialize the Parser object.

        Args:
            tf_api: an active instance of Text-Fabric with
                corpus source_data pre-loaded
        """

        # set up TF methods
        self._tfapi = tf_api
        self._F = tf_api.F
        self._E = tf_api.E
        self._T = tf_api.T
        self._L = tf_api.L

        # precedence assignments for when multiple matches are
        # found higher value == the function will be matched first
        self._prec = {
        }

    def _getprec(self, name):
        """Get precedence for parse values."""
        return self._prec.get(name, 0)

    def _getparsers(self):
        """Retrieve defined parsers in this class.

        Parsers are identified by those methods which
        are named without a prefixed underscore.

        Returns:
            dict where key=parsing value and 
            the value is the evaluation function
            which takes a single argument, a TF node,
            and returns True if a match is found 
        """

        # retrieve defined parsers
        parsenames = [
            m for m in dir(self) 
                if not m.startswith('_')
        ]

        # sort by precedence
        parsenames = sorted(
            parsenames,
            key=self._getprec,
            reverse=True,
        )

        # remap method names to actual functions for calls
        parsers = {
            p:getattr(self, p) for p in parsenames
        }
        
        return parsers

    def _parse(self, node):
        """Parse a supplied BHSA word node."""

        # get parsers sorted by precedence
        parsers = self._getparsers()

        # return first match in order of precedence
        for value, parser in parsers.items():
            if parser(node):
                return value

    def _getP(self, node, context='phrase_atom'):
        """Get Positions object for a TF node."""
        return PositionsTF(node, context, self._tfapi).get
