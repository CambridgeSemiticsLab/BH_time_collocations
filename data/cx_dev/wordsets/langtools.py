"""
Language Tools with Text-Fabric

This module contains classes to be used
with a Text-Fabric corpus instance.
"""

class Dummy:
    """A place-holder Class for use when supplied node is None"""
    def __init__(self, *args, **kwargs):
        pass
    def __bool__(self):
        return False
    def __getattr__(self, name):
        return lambda *args, **kwargs: None

class Positions:
    """Access positions around a node in a context.
    
    For a given origin node, provides data on another node
    that is (+/-)N positions away in a context.
    """
    
    def __init__(self, n, context, tf, order='slot'):
        """Prepare context and positions for a supplied TF node.
        
        Arguments:
            n: an integer that is a Text-Fabric node that serves
                as the origin node.
            context: a string that specifies a context in which
                to search for a position relative to origin node.
            tf: an instance of Text-Fabric with a loaded corpus.
            order: The method of order to use for the search.
                Options are "slot" or "node."
        """
        
        self.tf = tf.api
        self.n = n
        self.method = order
        self.thisotype = tf.api.F.otype.v(n)
        context = tf.api.L.u(n, context)[0]
        self.positions = tf.api.L.d(context, self.thisotype)     
        if order == 'node':
            self.originindex = self.positions.index(n)
    
    def nodepos(self, position):
        """Get position using node order.
        
        !CAUTION!
            This method should only be used with
            linguistic units known to be non-overlapping. 
            A "slot" object is a good example for which this 
            method might be used. 
            
            For example, given a phrase with another embedded phrase:
            
                > [1, 2, 3, [4, 5], 6]
                
            this method will indicate that slots 3 and 6 are adjacent
            with respect to the context. This is OK because we know 
            3 and 6 do not embed one another.
            By contrast, TF locality methods will mark 3 and 4 as adjacent.
            
        Arguments: 
            position: integer that is the position to find
                from the origin node
        """
        # use index in positions to get adjacent node
        # return None when exceeding bounds of context
        pos_index = self.originindex + position
        pos_index = pos_index if (pos_index > -1) else None
        try:
            return self.positions[pos_index]
        except (IndexError, TypeError):
            return None
    
    def slotpos(self, position):
        """Get position using slot order.
        
        This method is ideal for nodes that may overlap.
        Uses TF locality methods. These methods use
        pre-calculated tables to define whether two
        nodes are adjacent or not. The TF definition 
        of order is as follows:
            > if two objects are intersecting, 
            > but none embeds the other, 
            > the one with the smallest slot that 
            > does not occur in the other, comes first.
            > (https://annotation.github.io/text-fabric/Api/Nodes/)
            
        For example, given a phrase with another embedded phrase:
            
            > [1, 2, 3, [4, 5], 6]
        
        this method will indicate that slots 3 and 4 are adjacent.
        
        Arguments: 
            position: integer that is the position to find
                from the origin node
        """
        
        L = self.tf.L

        # determine which TF method to use
        if position < 0:
            move = L.p
        else:
            move = L.n
            
        # iterate the number of steps indicated by position
        # call move method for each step and set result to next position
        get_pos = self.n
        for count in range(0, abs(position)):
            get_pos = next(iter(move(get_pos, self.thisotype)), 0)
            
        if get_pos in self.positions:
            return get_pos
        else:
            return None
    
    def get(self, position, *features):
        """Get data on node (+/-)N positions away. 
        
        Arguments:
            position: a positive or negative integer that 
                tells how far away the target node is from source.
            features: a feature string or set to return based
                from the target node. If features not specified,
                will return the node itself.
        """
            
        # get next position based on method
        if self.method == 'slot':
            get_pos = self.slotpos(position)
        elif self.method == 'node':
            get_pos = self.nodepos(position)
    
        # return requested data
        if get_pos:
            Fs = self.tf.Fs
            if not features: 
                return get_pos
            elif len(features) == 1:
                return Fs(features[0]).v(get_pos)
            elif len(features) > 1:
                return set(Fs(feat).v(get_pos) for feat in features)
            
        # return empty data
        elif get_pos not in self.positions:
            if not features:
                return None
            elif len(features) == 1:
                return ''
            elif len(features) > 1:
                return set()

class Walker:    
    """Prepares paths from a source TF node to a target node.
    
    Supplies methods to walk forward or backward in a context until 
    encountering a TF node that meets a set of conds.
    
    Methods:
        ahead: Walk ahead from start node to target.
        back: Walk back from start node to target.
        firstresult: Return the first node in a path that
            returns True for a supplied function
    """
    
    def __init__(self, n, context, tf=None):
        """Initialize paths for a node.
        Arguments:
            n: Text-Fabric corpus node
            context: otype string of the supplied node's context to lookup
            tf: Running instance of Text-Fabric corpus
        """
        tf = tf.api      
        thisotype = tf.F.otype.v(n)
        context = tf.L.u(n, context)[0]
        self.positions = list(tf.L.d(context, thisotype))
        self.index = self.positions.index(n)

    def ahead(self, val_funct, **kwargs):
        """Walk ahead to node.
    
        Returns:
            Integer which corresponds to a Text-Fabric node or 
            output of function if output=True.
            
        Args:
            val_funct: a function that accepts a node argument
                and returns Boolean. This determines which word
                to return.          
                
        *Kwargs:
            stop: a function that accepts a node argument and
                returns Boolean. Determines whether to interrupt 
                the walk and return None.
            go: opposite of stop, a function that accepts a node
                argument and returns Boolean. Determines whether
                to keep going in a walk.
            output: return output of the val_funct instead of the
                node itself.
        """
        path = self.positions[self.index+1:]
        return self.firstresult(path, val_funct, **kwargs)
            
    def back(self, val_funct, **kwargs):
        """Walk back to node.
        
        Returns:
            Integer which corresponds to a Text-Fabric node or 
            output of function if output=True.
            
        Args:
            val_funct: a function that accepts a node argument
                and returns Boolean. This determines which word
                to return.          
                
        *Kwargs:
            stop: a function that accepts a node argument and
                returns Boolean. Determines whether to interrupt 
                the walk and return None.
            go: opposite of stop, a function that accepts a node
                argument and returns Boolean. Determines whether
                to keep going in a walk.
            output: return output of the val_funct instead of the
                node itself.
        """
        path = self.positions[:self.index]
        path.reverse()
        return self.firstresult(path, val_funct, **kwargs)
        
    def firstresult(self, path, val_funct, **kwargs):
        """Return data on node in a loop.
        
        Args:
            path: a set of nodes to traverse.
            val_funct: a function that accepts a TF node argument
                and returns Boolean. Determines which word to return
                in the walk.
                
        *Kwargs:
            stop: a function that accepts a node argument and
                returns Boolean. Determines whether to interrupt 
                the walk and return None.
            go: opposite of stop, a function that accepts a node
                argument and returns Boolean. Determines whether
                to keep going in a walk.
            output: return output of the val_funct instead of the
                node itself.
        """
        stop = kwargs.get('stop') or (lambda n: False)
        go = kwargs.get('go') or (lambda n: True)
        for node in path:
            # do matches
            test = val_funct(node)
            if test:
                if not kwargs.get('output', False):
                    return node
                else:
                    return test
            # do interrupts on go
            elif not go(node):
                break
            # do interrupts on stop
            elif stop(node):
                break
   