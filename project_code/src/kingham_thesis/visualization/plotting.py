import seaborn as sns
from bidi.algorithm import get_display

def format_rtl(string):
    """Format RTL Hebrew text for visualization."""
    return (
        get_display(string)
            .replace('\u05c1', '') # rm shin dots 
            .replace('\u05c2', '') 
    )    

def heatmap(data, center=0, **kwargs):
    """Draw seaborne heatmap with custom settings"""
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    hmkwargs = dict(
        center=center,
        cmap=cmap,
        square=True,
        linewidth=0.5   
    )
    hmkwargs.update(**kwargs)
    sns.heatmap(
        data,
        **hmkwargs,
    )

