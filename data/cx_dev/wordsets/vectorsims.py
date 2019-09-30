# import fasttext
# from scipy.spatial.distance import cosine
# sem_model = os.path.expanduser('~/github/codykingham/bhsa_vectors/model.bin') # semantic data


        
#         # !! Add cosine methods for pairwise comparisons
#         self.sim = fasttext.load_model(sem_model)
        
        
#     # !! For semantic vector model evaluations, must build in somehow
#     good_chars = set(char for word in F.otype.s('word') for char in F.voc_lex_utf8.v(word))
#     good_chars.remove('־') # remove maqqeph
#     def rem_accent(word):
#         # Remove accents from words
#         return ''.join(c for c in word if c in good_chars)