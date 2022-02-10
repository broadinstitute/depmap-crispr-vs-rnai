
import argparse
import sys
sys.path.append('/tmp/pipeline/src/dependency_metrics/packages/Python')
import pandas as pd
import numpy as np
import inference
import utilities as ut

def getParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-effect',
        type=str,
        default='none',
        help='File containing gene scores')
    parser.add_argument('--tpm',
        type=str,
        default=None,
        help='File containing RNAseq TPM expression data')
    parser.add_argument('--common-essential',
        type=str,
        default=None,
        help='Table with genes as rownames and boolean Common_Essential column')
    parser.add_argument('--pos-cntrl',
        type=str,
        default='none',
        help='Table with single column of gene names')
    parser.add_argument('--neg-cntrl',
        type=str,
        help='Table with single column of gene names')
    parser.add_argument('--pr-outfile',
        type=str,
        help='File for writing probabilities')
    parser.add_argument('--fdr-outfile',
        type=str,
        help='File for writing FDRs')
    return parser


def infer(data, initial_lambdas=[.85, .15], smoothing=.05, p_smoothing=.15):
    '''
    For each line, assign each gene score in that line a probability that the gene is a real dependency
    Parameters:
        data: melted gene frame with columns
             'lines': cell line name
              'genes': gene name
              'score': CERES score (scaled)
              'essential': Boolean whether the gene is identified essential (all lines)
              'unexpressed': Boolean whether the gene is expressed in that cell line
              'nonessential': Boolean whether the gene is identified nonessential (all lines)
        initial_lambdas: list of two float summing to 1;
                         best initial guess for how much of the distribution is null vs real dependency
        smoothing, p_smoothing: kernel parameters for distribution estimation, see cds.inference.probability_2class
    Returns:
        gene_dependency: matrix (DataFrame) of probabilities of each gene being a dependency in each cell line
        gene_fdr: matrix (DataFrame) of FDR values (see API)
    '''
    data['p_dependent'] = np.nan
    data['FDR'] = np.nan
    for line in data.lines.unique():
        inference.mixture_model.assign_probabilities(data, line, initial_lambdas, fit_gaussian=False,
                                         smoothing=smoothing, p_smoothing=p_smoothing)
    data.sort_values(by='p_dependent', ascending=False, inplace=True)
    for line in pd.unique(data.lines):
        vals = np.cumsum(
            1-data[data.lines == line].p_dependent
        )
        data.loc[vals.index, 'FDR'] = vals*1.0/np.arange(1, len(vals)+1)

    assert data.p_dependent[data.unexpressed].mean() < .2, 'mean unexpressed dependency is %f' %(
        data.p_dependent[data.unexpressed].mean())
    assert data.p_dependent[data.nonessential].mean() < .2, 'mean nonessential dependency is %f' %(
        data.p_dependent[data.nonessential].mean())
    # assert data.p_dependent[data.essential].mean() > .6, 'mean essential dependency is %f' %(
    #     data.p_dependent[data.essential].mean())

    gene_dependency = data.pivot(index='lines', columns='genes', values='p_dependent')
    gene_dependency.index.name = 'line'
    gene_fdr = data.pivot(index='lines', columns='genes', values='FDR')
    gene_fdr.index.name = 'line'
    return gene_dependency, gene_fdr

def melt_and_annotate(gene_effect, expression, essentials, negative_controls, expression_cutoff=.2):
    '''
    produces a melted DataFrame with annotations necessary for assigning dependency probabilities
    Parameters:
        gene_effect: matrix with genes as columns
        expression: matrix of TPM expression values with genes as columns
    Returns:
        melted DataFrame
    '''

    gene_effect = gene_effect.copy()
    #long form of gene_effect
    melted = ut.index_melt(gene_effect) 
    try:
        melted.columns = ['lines', 'genes', 'score']
        #add column with boolean if gene is in essentials list
        melted['essential'] = melted.genes.isin(essentials) 
        #add column with boolean if gene is in negative_controls list
        melted['nonessential'] = melted.genes.isin(negative_controls)
        assert sum(melted.essential) > 0 and sum(melted.nonessential > 0), "no genes found that match the positive or negative control sets"
    except Exception as e:
        print(melted[:10])
        print(essentials[:10])
        print(negative_controls[:10])
        raise e

    if not expression is None:
        #Melt the expression matrix to match the gene_effect matrix
        expression = expression.reindex(gene_effect.index)
        expression = expression.reindex(columns=gene_effect.columns)
        melt_exp = ut.index_melt(expression)
        #Append column to melted gene_effect to indicate if gene is not expressed
        melted['unexpressed'] = melt_exp.value < expression_cutoff
    else:
        melted['unexpressed'] = melted['nonessential']

    #Assert that positive and negative control groups have reasonable means
    mu = melted[melted.unexpressed].score.mean()
    if not np.abs(mu) < .1:
        raise RuntimeError("Unexpressed mean displaced from 0 (%f)" %(
        np.abs(mu) ))
    mu2 = melted[~melted.unexpressed].score.mean()
    if not mu2 < mu:
        raise RuntimeError("Unexpressed mean left of expressed mean (%f vs %f)" %(
        mu, mu2 ))
    mu3 = melted[melted.essential].score.mean()
    #if not mu3 < -.6:
    #    raise RuntimeError('Median of essentials shifted to %f' %mu3)
    if not mu3 < mu:
       raise RuntimeError('Median of essentials shifted to %f' %mu3)
    mu4 = melted[melted.nonessential].score.median()
    if abs(mu4) > .02:
        raise RuntimeError('Median of nonessentials shifted to %f' %mu4)

    melted.unexpressed.fillna(False, inplace=True)
    melted = melted[~melted.score.isnull()]
    return melted

if __name__ == '__main__':
    args = getParser().parse_args()

    #Load gene effect from csv where first column is row names
    gene_effect = pd.read_csv(args.gene_effect,index_col=0)

    #Load expression data and filter by genes included in gene effect
    if not args.tpm is None:
        exp = pd.read_csv(args.tpm,index_col=0)
        exp = exp.loc[exp.index.intersection(gene_effect.index),exp.columns.intersection(gene_effect.columns)]
    else:
        exp = None

    #common essentials and unexpressed have priority over pos and neg control if files provided
    if not args.common_essential is None:
        essentials = pd.read_csv(args.common_essential,index_col=0)
        essentials = essentials.index[essentials.Common_Essential]
    else:
        essentials = pd.read_csv(args.pos_cntrl,index_col=0)
        essentials = essentials.index.tolist()
    
    negative_controls = pd.read_csv(args.neg_cntrl,index_col=0)
    negative_controls = negative_controls.index.tolist()

    #Annotate the dataframe with which genes are non-expressed in TPM data
    melted  = melt_and_annotate(gene_effect, exp, essentials, negative_controls, expression_cutoff=.2)

    #Run the inference method
    gene_dependency, gene_fdr = infer(melted, initial_lambdas=[.85, .15], smoothing=.05, p_smoothing=.15)

    gene_dependency.index.name = 'Row.name'
    gene_dependency.to_csv(args.pr_outfile,index=True)

    gene_fdr.index.name = 'Row.name'
    gene_fdr.to_csv(args.fdr_outfile,index=True)










