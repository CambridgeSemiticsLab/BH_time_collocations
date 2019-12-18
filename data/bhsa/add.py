"""
In this module, new nodes and features
will be added to the BHSA dataset
"""

from tf.fabric import Fabric
from tf.compose import modify
from timephrase import chunk_time

def to_graph(locations, base_metadata):
    """
    Add new nodes to the custom TF dataset.
    """

    # setup data and methods
    bhsa, output = locations['bhsa'], locations['custom']
    TF = Fabric(locations=[bhsa,output], silent='deep') 
    api = TF.load('function language')
    F, L, N = api.F, api.L, api.N

    # set up existing BHSA graph data 
    nodeFeatures =  {
        'otype': {
            n:F.otype.v(n) for n in N() 
        },
    }
    edgeFeatures =  {
        'oslots': {
            n:L.d(n,'word') for n in N()
                if F.otype.v(n) != 'word'
        },
    }
    metaData = {
        '': base_metadata,
        'oslots': {'valueType':'int','edgeValues':False},
        'otype': {'valueType':'str'},
    }

    # populate new node data here
    addTypes = {}

    # populate new node data
    chunk_time(addTypes, api)

    # add new objects and features 
    maxNode = max(nodeFeatures['otype'])
    nodeMap = {}
    for type, data in addTypes.items():
        for node, slots in data['nodeSlots'].items():
            maxNode += 1
            nodeFeatures['otype'][maxNode] = type
            edgeFeatures['oslots'][maxNode] = slots
            nodeMap[node] = maxNode

        # add aditional features if present
        for feat, data in addTypes[type].get('nodeFeatures',{}).items():
            nodeFeatures[feat] = {}
            for node, value in data.items(): 
                nodeFeatures[feat][nodeMap[node]] = value

        # option: insert edge feature collection here

        # populate meta data
        metaData.update(addTypes[type].get('metaData',{}))

    # export data
    TFs = Fabric(locations=output, silent=True)
    TFs.save(
        nodeFeatures=nodeFeatures,
        edgeFeatures=edgeFeatures,
        metaData = metaData,
        silent=True,
    )    
