import json
from tf.fabric import Fabric
from pos_parser import posParser

# initialize Text-Fabric
tf = Fabric(locations=snakemake.input.bhsadata)
tf_api = tf.load('''
    pdp sp ls lex st 
    prs function vt
''')

F = tf_api.F

# assign parts of speech
# if posParser does not return a value, we 
# assign it the default BHSA tag in uppercase
parser = posParser(tf_api)

slot2pos = {}
for slot in F.otype.s('word'):
    pos = parser._parse(slot)
    pos = pos or F.pdp.v(slot).upper()
    slot2pos[slot] = pos

uniquepos = sorted(set(slot2pos.values()))

# export to JSON
with open(snakemake.output.slot2pos, 'w') as outfile:
    json.dump(slot2pos, outfile, indent=1)

# export unique values, which are needed to configure the phrase parser
with open(snakemake.output.uniquepos, 'w') as outfile:
    outfile.write(',\n'.join(uniquepos))
