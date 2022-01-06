from __future__ import print_function
import numpy as np
import pandas as pd
#from taigapy import TaigaClient

def covar(A, B=None, axis=0):
    #adapted from https://stackoverflow.com/a/30143754/8691282
    if B is None:
        B = A
    if axis == 0:
        A = A.T
        B = B.T
    elif axis != 1:
        raise ValueError('axis must be 0 or 1, not %r' %axis)
    A_mA = A.subtract(A.mean(axis=1), axis=0)
    A_mA[A_mA.isnull()] = 0
    B_mB = B.subtract(B.mean(axis=1), axis=0)
    B_mB[B_mB.isnull()] = 0

    # Sum of squares across rows
    ssA = (A_mA**2).sum(axis=1).as_matrix()
    ssB = (B_mB**2).sum(axis=1).as_matrix()

    # Finally get corr coeff
    out = pd.DataFrame(np.dot(A_mA.as_matrix(),B_mB.T.as_matrix())/np.sqrt(np.dot(ssA[:,None],ssB[None])),
                        index=A.index, columns=B.index)
    return out


def roll(x, dx):
    '''allow rolling numpy arrays by floats rather than integer indices'''
    if dx // 1 and dx % 1:
        return roll(roll(x, dx // 1), dx % 1)
    elif dx // 1:
        return np.roll(x, int(dx))
    elif dx % 1:
        out = (1 - abs(dx % 1)) * x + abs(dx % 1) * np.roll(x, int(np.sign(dx)))
        return out
    else:
        return x


def differential(x, step, dx=1):
    '''Find the derivative of an array using periodic boundary conditions'''
    return (roll(x, -.5 * dx) - roll(x, .5 * dx)) * 1.0 / (dx * step)


def smooth(x, sigma):
    '''smooth with a gaussian kernel'''


    kernel = np.exp(-.5 * np.arange(int(-4 * sigma), int(4 * sigma + 1), 1) ** 2 / sigma ** 2)
    kernel = kernel / sum(kernel)
    return np.convolve(x, kernel, 'same')


def get_mode(vec, nbins=None):
    '''
    heuristic method for getting mode of a set of points
    '''
    vec = np.array(vec)
    vec = vec[~np.isnan(vec)]
    if nbins is None:
        nbins = int(len(vec)/25)
    heights, bounds = np.histogram(vec, bins=nbins)
    max_ind = np.argmax(heights)
    if max_ind < 2 or max_ind > len(bounds) - 4:
        lower = max_ind
        upper = max_ind + 1
    else:
        lower = max_ind - 2
        upper = max_ind + 3
    try:
        out = np.median(vec[(vec >= bounds[lower]) & (vec <= bounds[upper])])
    except IndexError:
        raise RuntimeError('%r\n%r\n%i %i %i %i' %(heights, bounds, max_ind, lower, upper, len(bounds)))
    return out


def get_binned(data, low=None, high=None, bin_width=None, bins=None):
    '''
    Get a pandas dataframe of the normalized histogram of the data
    args:
        data: 1D array or pandas Series
        low, high: left- and right-most boundaries of the bins
        bin_width: if specified, bins will have this width (float)
        bins: if specified and bin_width is not, bins will be sized to have this many total bins
    returns:
        pandas DataFrame with x holding bin centers and y the normalized counts in each bin
    '''
    if low is None:
        low = min(data)
    if high is None:
        high = max(data)
    if bin_width:
        breaks = np.arange(low, high + bin_width, bin_width)
    elif bins:
        breaks = np.linspace(low, high, bins + 1)
    else:
        raise ValueError("One of bin_width or bins must be supplied")
    x = breaks[:-1] + .5 * bin_width
    bins = pd.DataFrame({'x': x})
    bins['y'] = np.histogram(data, breaks, density=True)[0]
    return bins


# pulling data

def get_crispr(name='avana-1.0'):
    ''' get the CERES matrix'''
    c = TaigaClient()
    solid_data = c.get(name=name)
    solid_data.rename(columns={solid_data.columns[0]: 'genes'}, inplace=True)
    return solid_data

def index_melt(df, **kwargs):
    name = df.index.name
    if not name:
        name = 'index'
    return pd.melt(df.reset_index(level=(df.index.nlevels-1)), id_vars=name, **kwargs)

def load_features(**kwargs):
    '''safer loading of gene/cell_line data that resolves common discrepancies. Kwargs are
    passed to TaigaClient.get'''
    c = TaigaClient()
    try:
        df = c.get(**kwargs)
    except TypeError:
        raise ValueError('no data found for %r' %kwargs)

    #check if cell lines are columns - this heuristic has been safe so far
    flipped = len(df.columns[3].split('_')) > 1
    #If the first column is not strictly numeric, we'll assume its an index
    try:
        _ = df[df.columns[0]].astype(np.float)
        column_index = False
    except ValueError:
        column_index = True
    if column_index:
        df.set_index(df.columns[0], inplace=True)
    if flipped:
        df = df.T
    #check if genes include Broad ID
    if len(df.columns[3].split('(')) > 1:
        df.columns = [s.strip().split(' ')[0] for s in df.columns]
    if any(df.index.duplicated()):
        raise RuntimeError('index duplicates: %r' %df.index[df.index.duplicated()])
    if any(df.columns.duplicated()):
        raise RuntimeError('column duplicates: %r' %df.columns[df.columns.duplicated()])
    return df

def get_melted(essentials=[], nonessentials=[], name='avana-564b', version=None):
    '''
    get a collection of relevant cell line X gene data in a melted form
    Parameters (optional):
        essentials: list of genes considered panessential
        nonessentals: list of genes considered nonessential
        name: taiga name of dataset for dependency data, defaults to latest
        version: taiga version
    '''
    print('loading scores')
    scores = load_features(name=name, version=version, file='gene_effect')
    print('loading probabilities')
    probabilities = load_features(name=name, version=version, file='gene_dependency')
    probabilities = probabilities.loc[scores.index, scores.columns]
    print('loading FDR')
    fdr = load_features(name=name, version=version, file='gene_fdr')
    fdr = fdr.loc[scores.index, scores.columns]
    print('loading expression')
    expressions = load_features(name=
                        "ccle-rnaseq-expression-genes")
    expressions = expressions.loc[scores.index, scores.columns]
    print('melting')
    melted = index_melt(scores)
    melted.columns = ['lines', 'genes', 'score']
    melted['p_dependent'] = index_melt(probabilities).value
    melted['expression'] = index_melt(expressions).value
    melted['unexpressed'] = melted.expression < -2.9
    melted['fdr'] = index_melt(fdr).value
    melted.unexpressed[pd.isnull(melted.unexpressed)] = False
    melted['essential'] = melted.genes.isin(essentials)
    melted['nonessential'] = melted.genes.isin(nonessentials)
    return melted

# GSEA file format savers
def to_gsea_cls(series, filename):
    classes = series.unique()
    class_map = pd.DataFrame({'class_name': classes, 'index': list(range(len(classes)))}).set_index('class_name')
    with open(filename, 'w+') as f:
        f.write('%i %i 1\n' %(len(series), len(classes)))
        f.write('#')
        for c in classes:
            f.write(' %s' %c)
        f.write('\n')
        for c in series.values:
            f.write('%i ' %class_map.loc[c])
            
def to_gsea_txt(dataframe, filename):
    out = dataframe.copy()
    out['DESCRIPTION'] = 'nan'
    out = out[['DESCRIPTION'] + list(dataframe.columns)]
    out.index = [s.split(' ')[0] for s in out.index]
    out.index.name = 'NAME'
    out.to_csv(filename, sep='\t')


# validation functions

def shuffle(df):


    '''shuffle each column of a dataframe independently'''
    indices = np.arange(len(df), dtype=np.int)
    shuffled = np.copy(df.as_matrix())
    for i in range(df.shape[1]):
        np.random.shuffle(indices)
        shuffled[:, i] = shuffled[indices, i]
    return pd.DataFrame(shuffled, index=df.index, columns=df.columns)


