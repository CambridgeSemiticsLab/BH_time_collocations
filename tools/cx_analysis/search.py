import random
import itertools
import networkx as nx
from IPython.display import display, HTML
from datetime import datetime
from pprint import pprint
from .cx import Construction

class SearchCX:
    """Methods for visualizing construction objects with TF"""
    
    def __init__(self, tf_app):
        """Set up TF methods for class-wide use"""
        api = tf_app.api
        self.F,self.T,self.L = api.F, api.T, api.L
        self.app = tf_app
        
    def pretty(self, obj, condense='phrase', **kwargs):
        """Show a linguistic object that is not native to TF app."""
        L = self.L
        A = self.app
        index = kwargs.get('index')
        kwargs = {k:v for k,v in kwargs.items() if k not in {'index'}}
        show = L.d(obj, condense) if index is None else (L.d(obj, condense)[index],)
        print(show, not index, index)
        A.prettyTuple(show, seq=kwargs.get('seq', obj), **kwargs)

    def prettyconds(self, cx):
        '''
        Iterate through an explain dict for a rela
        and print out all of checked conditions.
        '''
        cx_tree = [
            n for n in nx.bfs_tree(cx.graph, cx)
                if type(n) == Construction
        ]

        for node in cx_tree:
            print(f'-- {node} --')
            for case in node.cases:
                print(f'pattern: {case.get("pattern", case["name"])}')
                for cond, value in case['conds'].items():
                    print('{:<30} {:>30}'.format(cond, str(value)))
                print()

    def showcx(self, cx, **kwargs):
        """Display a construction object with TF.

        Calls TF.show() with HTML highlights for 
        words/stretch of words that serve a role
        within the construction. 
        """

        L = self.L
        A = self.app
        
        # get slots for display
        refslots = cx.slots if cx.slots else cx.element.slots
        showcontext = tuple(set(L.u(s, 'phrase')[0] for s in refslots))
        timephrase = L.u(list(refslots)[0], 'timephrase')[0]        

        if not cx:
            print('NO MATCHES')
            print('-'*20)
            A.prettyTuple(
                showcontext, extraFeatures='sp st', 
                withNodes=True, seq=f'{timephrase} -> {cx}'
            )
            if kwargs.get('conds'):
                self.prettyconds(cx)
            return None

        colors = itertools.cycle([
            '#96ceb4', '#ffeead', '#ffcc5c', '#ff6f69',
            '#bccad6', '#8d9db6', '#667292', '#f1e3dd',
        ])
        highlights = {}
        role2color = {}

        for node in cx.graph.adj[cx]:
            role = cx.graph[cx][node]['role']
            slots = cx.getslots(node)
            color = next(colors)
            role2color[role] = color
            for slot in slots:
                highlights[slot] = color

        A.prettyTuple(
            showcontext, 
            extraFeatures=kwargs.get('extraFeatures', 'sp st lex'), 
            withNodes=True, 
            seq=f'{timephrase} -> {cx}', 
            highlights=highlights
        )
        # reveal color meanings
        for role,color in role2color.items():
            colmean = '<div style="background: {}; text-align: center">{}</div>'.format(color, role)
            display(HTML(colmean))

        pprint(cx.unfoldroles(), indent=4)
        print()
        if kwargs.get('conds'):
            self.prettyconds(cx)
        display(HTML('<hr>'))

    def search(self, elements, cxtest, pattern='', 
                    show=None, end=None, shuffle=True,
                    updatei=1000, select=None, **kwargs):
        """Search phrases for a specified relation"""

        start = datetime.now()
        print('beginning search')

        # random shuffle to get good diversity of examples
        if shuffle:
            random.shuffle(elements)
        matches = []

        # iterate and find matches on words
        for i,el in enumerate(elements):

            # update every 5000 iterations
            if i%updatei == 0:
                print(f'\t{len(matches)} found ({i}/{len(elements)})')

            # run test for construction
            test = cxtest(el)

            # save results
            if test:
                if pattern:
                    if test.pattern == pattern:
                        matches.append(test)
                else:
                    matches.append(test)

            # stop at end
            if end and len(matches) == end:
                break

        # display
        print('done at', datetime.now() - start)
        print(len(matches), 'matches found...')
        if show:
            print(f'showing top {show}')

        # option for filtering results
        if select:
            matches = [m for m in matches if select(m)]
            print(f'\tresults filtered to {len(matches)}')

        for match in matches[:show]:
            self.showcx(match, **kwargs)

        return matches
            
# NB: For the future. Here is a template to plot 
# a network graph using networkx.

# graph = GIVE GRAPH HERE

# plt.figure(figsize=(10,5))
# pos = nx.drawing.spectral_layout(graph)
# nx.draw_networkx(graph, pos)

# edge_labels = {
#     (n1,n2):graph[n1][n2]['role']
#         for n1,n2 in graph.edges
# }
    
# nx.draw_networkx_edge_labels(graph, pos, font_size=10, edge_labels=edge_labels)
# plt.show()