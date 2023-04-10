import json
from tf.fabric import Fabric

def patch_function(paths):
    """Remap phrase functions in BHSA.

    Patches are applied to phrase functions deemed problematic.
    """

    # setup data and methods
    TF = Fabric(locations=paths['bhsadata'], silent='deep')
    api = TF.load('function')
    F = api.F

    # features to be modded
    functions = {
        str(n):F.function.v(n)
            for n in F.otype.s('phrase')
    }

    # remap functions
    with open(paths['editfuncts'], 'r') as infile:
        editfuncts = json.load(infile) 

    for node, edit_data in editfuncts.items():
        functions[node] = edit_data['function']

    with open(paths['patched'], 'w') as outfile:
        json.dump(functions, outfile, indent=1)
