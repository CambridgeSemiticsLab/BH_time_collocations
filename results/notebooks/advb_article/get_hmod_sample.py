from __main__ import *
hm_df = functs_df[~((functs_df.head_type == 'prep') & (functs_df.suffix))].copy()
