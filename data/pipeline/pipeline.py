
# coding: utf-8

# # Pipeline
# 
# This module describes in detail the data-production pipeline of the Time Collocations project. The root source of the data is the [BHSA](https://www.github.com/etcbc/bhsa) of the Eep Talstra Centre for Bible and Computer, Vrije Universiteit Amsterdam. The BHSA contains grammatical and syntactic annotations, including a feature called `Function`, which identifies adverbial time phrases.
# 
# This project makes four primary kinds of modifications to the root BHSA data:
# 
# 1. **Some phrases marked for time adverbial function are corrected and remapped to the appropriate function.** These corrections are based on examining the data manually and deciding whether the default BHSA label should be kept or modified. Those cases will be described explicitly below.
# 2. **Meaningful groups above or below a phrase or other formal object are "chunked" into new objects, called `chunk`.** The `chunk` object contains word groups such as quantifier number chains (below the phrase level) or time adverbial chunks (units above the phrase level). This latter case consists of certain situations where BHSA splits a time adverbial phrase into two separate phrases. Those two parts are recombined in `chunks` for further analysis.
# 3. But will be number 4. **New features on phrases, chunks, and constructions.** These are various features stored on the objects noted above.
# 
# 4. **~TO BE DONE, will be 3.~** **Construction objects are generated**. This will be the final step in the pipeline, which is based on the statistical and functional analysis of time adverbials. A series of new objects will be built using the `chunk` objects. They will contain labels for semantic roles and semantic function. These labels are to be constructed using inductive, data-driven analysis combined with insights from Construction Grammar.
# 
# **The end result is a custom database, which consists of original BHSA data, modified BHSA phrase function data, and new objects with their associated features**.
# 
# The whole pipeline is represented below in diagram form. 
# 
# <table>
# <img src="../../docs/images/pipeline_diagram.png" width="30%" height="30%" align="middle">
# </table>

# # 0. Modules, Paths, and Classes
# 
# Input and output files are outlined below. Both directories consist (or will) of .tf resources. The input files are the base BHSA dataset while the output will be the customized, project dataset. Another dataset, `heads`, which is crucial for this work, is identified below as well. See the documentation for heads [here](https://nbviewer.jupyter.org/github/ETCBC/heads/blob/master/phrase_heads.ipynb).

# In[1]:


# import generic packages and modules
import os, sys, collections, glob
from tf.fabric import Fabric
from tf.app import use

# import pipeline classes
from remapfunctions import RemapFunctions # remaps phrase functions
from chunking import Chunker # chunks meaningful groups
from enhance import Enhance # adds helper data to chunks
from function_association import FunctAssoc # calcs function associations

# function associations take a long time to run
# only run them if explicitly asked
run_associations = '-full' in sys.argv or False # config to True in NB to run it

# determine whether script is being run from ipython or as simple .py
# this notebook code is converted with nbconvert to a .py script and has
# additional visualizations that do not need to be run in the .py version.
try:
    is_nb = __IPYTHON__
except: 
    is_nb = False
    
# load Text-Fabric instance if running inside NB
# For visualizing data that will be exported
if is_nb:
    bhsa = use('bhsa', silent=True)
    F, T, L = bhsa.api.F, bhsa.api.T, bhsa.api.L
    
# configure input TF data and output dir for new TF data
home = os.path.expanduser('~/')
bhsa_dir = os.path.join(home, 'text-fabric-data/etcbc/bhsa/tf/c') # input
heads_dir = os.path.join(home, 'github/etcbc/heads/tf/c') # heads data
output_dir = '../tf' # output; all new data goes here
locations = {'bhsa': bhsa_dir,
             'heads': heads_dir,
             'output': output_dir}

# The metadata below is assigned to all new features.
base_metadata = {'source': 'https://github.com/etcbc/bhsa',
                 'data_version': 'BHSA version c',
                 'origin': 'Made by the ETCBC of the Vrije Universiteit Amsterdam; edited by Cody Kingham'}


# ### Purge Old .tf Files
# 
# Remove all .tf files in the output directory. This is necessary since all output data goes to the same place, and subsequent classes depend on previous data. Without the purge, subsequent runs will add new objects on top of previous ones!
# 
# The `funct_assoc.tf` and `top_assoc.tf` features rely only on the phrase function feature and need not be purged unless new phrase function edits are added (cf. #1). Keeping them prevents the significant runtime needed to calculate them.

# In[2]:


# keep association data by default
# these can be written over without any problems
keep = [os.path.join(output_dir, file) for file in ('funct_assoc.tf', 'top_assoc.tf')] 

print('Purging old data...')
old_tf = glob.glob(os.path.join(output_dir, '*.tf'))
for file in old_tf:
    if file not in keep:
        print(f'\tpurging {file}')
        os.remove(file)
print('DONE')


# # 1. Phrase Function Edits

# In[3]:


remap_functs = RemapFunctions(locations, base_metadata)


# The following functions will be changed...

# In[4]:


if is_nb:
    for phrase, newfunct in remap_functs.newfunctions.items():
        print(f'{phrase} node with function {F.function.v(phrase)} will be remapped to {newfunct}')
        bhsa.pretty(phrase)


# The change is executed below.

# In[5]:


remap_functs.execute()


# # 2.1 Generate Chunk Objects
# 
# The ETCBC data is not granular enough for many types of searches beneath the phrase level. For example, if there are coordinated noun phrases that all function as a single phrase, the individual phrases are not delineated. I need them spliced out so that I can track coordinated nouns in a phrase. Another case is with quantifiers, wherein the quantifier chains themselves are not in any way set apart from other items in the phrase. For these cases, I will make `chunk` objects—these are essentially phrase-like objects. Another problem is that some phrases are split into two, whereas elsewhere in the database the same phrase pattern is portrayed as a single phrase. This is fixed by creating a new `chunk` object.
# 
# `chunk` objects have two important features, called `label` and `role`. The `label` feature is essentially the name of the function. For instance, `timephrase` or `quant` (quantifier). The feature `role` is an edge feature, which maps the chunk's component parts to the new object. For example, in a `quant_NP`, there is an edge drawn from the noun to the `quant_NP` chunk; the edge has a value of "quantified" since the noun is the quantified item.

# In[6]:


chunker = Chunker(locations, base_metadata)


# ### [relevant, illustrative examples will be shown here]

# In[7]:


chunker.execute()


# # 3. Generate Helper Data
# 
# It is easier to use the resulting TF resource to build features of chunks than building those features from the raw data. In this section, enhancements are added to the chunk objects that were generated in section 2.

# In[8]:


enhance = Enhance(locations, base_metadata)


# ## 3.1 Embedding Data
# 
# Quantifier chunks often consist of smaller, component chunks. When clustering time adverbials, it is often not necessary to know that a quantifier contains two component parts. Rather, only the top level quantifier chunk is needed to indicate that there is quantification in the adverbial. The `embed` feature is a simple `true` or `false` tag which indicates whether a given chunk is contained within another chunk of the same kind. This allows quick and efficient selection of non-embedded chunks by `embed=false`. 

# In[9]:


enhance.embeddings()


# ## 3.2 Quantifier Time Roles
# 
# Quantifiers contained in a timephrase do not yet have a role mapping from the quantified noun to the time chunk. We add that below, creating an edge feature with a role of 'time' from a quantified noun to its embedding `timephrase` chunk.

# In[10]:


enhance.quanttimes()


# ## Export Enhancements
# 
# The enhancements are now exported to .tf resources.

# In[11]:


enhance.export()


# ## 4. Function Associations
# 
# Build association scores between head lexemes and their phrase functions. This data takes a significant amount of time to recalculate and only needs to be run if function data has changed.

# In[5]:


if run_associations:
    funct_assoc = FunctAssoc(locations, base_metadata)
    funct_assoc.execute()

