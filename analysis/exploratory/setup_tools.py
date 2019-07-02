import os
from importlib import util
import collections
import csv
import random
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(font_scale=1.5, style='whitegrid')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from skfuzzy.cluster import cmeans
from tf.fabric import Fabric
from tf.app import use

# import custom tools
os.sys.path.append('/Users/cody/github/CambridgeSemiticsLab/time_collocations/analysis/')
from tools.locations import data_locations 
from tools.significance import contingency_table, apply_fishers
from tools.pca import plot_PCA
from tools.helpers import convert2pandas, funct2function, formatPassages
from tools.visualize import reverse_hb, countBarplot
from tools.tokenizers import tokenize_surface
from tools.time import Time

TF = Fabric(locations=data_locations.values())
api = TF.load('''

vs vt pdp gloss lex language ps gn
rela typ number function prs
g_cons_utf8 lex_utf8 nu mother st uvf
g_word_utf8 trailer_utf8 voc_lex_utf8
head nhead obj_prep sem_set
ls top_assoc funct_assoc kind txt
label role
''')

A = use('bhsa', api=api, hoist=globals(), silent=True)

A.displaySetup(condenseType='clause', condensed=True, withNodes=True)

firstyear = '../../data/paper_data/firstyear2' # directory for first year review paper saves