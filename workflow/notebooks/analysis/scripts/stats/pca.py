'''
Module contains a useful method for plotting PCA data with tags annotated.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns

def apply_pca(df, sample_axis, feature_axis,
              scree=True, components=0, scale=True):
    """Apply PCA analysis to a supplied dataset.

    Args:
        df: dataframe with data
        sample_axis: 0 (rows) or 1 (columns), axis of samples
        feature_axis: 0 (rows) or 1 (columns), axis of features 
        scree: whether to display scree plot of explained
            variance
        components: number of components to calculate
        scale: whether to normalize the values

    Returns:
        (df_with_pca_values, df_with_feature_loadings)
    """

    # determine size/shape of data
    # NB from sklearn documentation: pca.fit expects shape of: 
    #     >>> "array-like, shape (n_samples, n_features)"
    # Thus we make sure to invert the table if necesssary
    if sample_axis == 1 and feature_axis == 0:
        df = df.T
    elif sample_axis != 0 and feature_axis != 1:
        raise Exception('Invalid axis! Should be 0 or 1')
    targets = df.index
    features = df.columns

    # preprocess
    if scale:
        scaler = StandardScaler()
        scaler.fit(df)
        data = scaler.transform(df)
    else:
        data = df.values

    # run the analysis
    n_components = components or features.shape[0]
    pca_analysis = PCA(n_components)
    pca_fit = pca_analysis.fit(data)
    components = pca_analysis.transform(data)
    pca_targets = pd.DataFrame(
        components,
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
        plt.figure(figsize=(6,5))
        explained_variance = pca_fit.explained_variance_ratio_[:n_components]
        sns.barplot(
            x = np.arange(n_components)+1,
            y = explained_variance,
            color='white',
            edgecolor='black',
        )
        plt.xlabel('Principle Component', size=12)
        plt.ylabel('Raio of Explained Variance', size=12)
        plt.title(
            f'Ratio of Explained Variance for Principle Components 1-{n_components}',
            size=12)
        plt.show()
        print('explained variance:')
        print(explained_variance)

    return (pca_targets, loadings, explained_variance)
