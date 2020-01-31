'''
Module contains a useful method for plotting PCA data with tags annotated.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns

def apply_pca(df, target_axis, feature_axis, scree=True, components=0):
    """Apply PCA analysis to a supplied dataset.

    Args:
        df: dataframe with data
        targets: 0 or 1, axis of target terms
        features: 0 or 1, axis of features for analysis
        scree: whether to display scree plot of explained
            variance
        components: number of components to calculate

    Returns:
        (df_with_pca_values, df_with_feature_loadings)
    """

    # determine size/shape of data
    # NB from sklearn documentation: pca.fit expects shape of: 
    #     >>> "array-like, shape (n_samples, n_features)"
    # Thus we make sure to invert the table if necesssary
    if target_axis == 1 and feature_axis == 0:
        df = df.T
    elif target_axis != 0 and feature_axis != 1:
        raise Exception('Invalid axis! Should be 0 or 1')
    targets = df.index
    features = df.columns

    # run the analysis
    n_components = components or features.shape[0]
    pca_analysis = PCA(n_components)
    pca_fit = pca_analysis.fit(df.values) 
    pca_values = pca_analysis.transform(df.values)
    
    pca_targets = pd.DataFrame(
        pca_values,
        columns = [f'PC{i+1}' for i in np.arange(n_components)],
        index = df.index
    )

    # compile loadings into DataFrame. There may be different 
    # conventions for loading calculations. I follow the 
    # definition and advice offered here:
    # https://stats.stackexchange.com/questions/143905/loadings-vs-eigenvectors-in-pca-when-to-use-one-or-another
    # that is: 
    # Loadings = Eigenvectors * sqrt(Eigenvalues)
    loadings = pca_fit.components_.T * np.sqrt(pca_fit.explained_variance_)
    loadings = pd.DataFrame(
        loadings.T, 
        index=np.arange(n_components)+1, 
        columns=df.columns
    )

    # plot optional scree-plot of explained variance for each component
    if scree:
        plt.figure(figsize=(8,6))
        sns.barplot(
            x = np.arange(n_components)+1,
            y = pca_fit.explained_variance_ratio_[:n_components],
            color='steelblue',
            edgecolor='black',
        )
        plt.xlabel('Principle Component', size=16)
        plt.ylabel('Raio of Explained Variance', size=16)
        plt.title(
            f'Ratio of Explained Variance for Principle Components 1-{n_components}',
            size=16)
        plt.show()

    return pca_targets, loadings

def plot_PCA(pca_nouns, 
             zoom=tuple(), 
             noun_xy_dict=False, 
             save='', 
             annotate=True, 
             title='', 
             components=tuple(),
             annoTags=[],
             anno_size='18'
            ):
    '''
    Plots a PCA noun space.
    Function is useful for presenting various zooms on the data.
    '''
    
    x, y = components
    
    # plot coordinates
    plt.figure(figsize=(12, 10))
    plt.scatter(x, y)

    if zoom:
        xmin, xmax, ymin, ymax = zoom
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
    
    if title:
        plt.title(title, size=18)
    plt.xlabel('PC1', size=18)
    plt.ylabel('PC2', size=18)
    plt.axhline(color='red', linestyle=':')
    plt.axvline(color='red', linestyle=':')
    
    # annotate points
    if annotate:
        noun_xy = {} # for noun_dict
        noun_lexs = annoTags
        
        for i, noun in enumerate(noun_lexs):
            noun_x, noun_y = x[i], y[i]
            noun_xy[annoTags[i]] = (noun_x, noun_y)
            if zoom: # to avoid annotating outside of field of view (makes plot small)
                if any([noun_x < xmin, noun_x > xmax, noun_y < ymin, noun_y > ymax]):                
                    continue # skip noun
            plt.annotate(noun, xy=(noun_x, noun_y), size=anno_size)
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    
    plt.show()
    
    if noun_xy_dict:
        return noun_xy
