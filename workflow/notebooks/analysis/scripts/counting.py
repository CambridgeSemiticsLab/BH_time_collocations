def pivot_ct(df, index, columns, sort=True, min_obs=0):
    """Fast count table."""
    ct = df.pivot_table(
        index=index,
        columns=columns,
        aggfunc='size',
        fill_value=0,
    )
    
    # sort
    if sort:
        ct = ct.loc[ct.sum(1).sort_values(ascending=False).index]
        ct = ct[ct.sum(0).sort_values(ascending=False).index]
        
    # filter
    ct = ct.loc[ct.sum(1) >= min_obs, :]
    ct = ct.loc[:, ct.sum(0) >= min_obs]
    
    return ct
