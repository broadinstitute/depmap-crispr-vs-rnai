from __future__ import print_function
import argparse
import sys
from time import time
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, ElasticNet
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.metrics import roc_auc_score, make_scorer, r2_score
from sklearn.feature_selection import SelectKBest, VarianceThreshold, mutual_info_regression, f_regression
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from scipy.stats import pearsonr
import json
import os.path
from pathlib import Path

import math
import random
from itertools import chain

import os
import gc




def getParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model-config',
        type=str,
        default='none',
        help='JSON file containing model definitions')
    parser.add_argument('--task-mode',
        type=str,
        default='regress',
        help='regress or classify')
    parser.add_argument('--targets',
        type=str,
        default='none',
        help='Matrix of dependency probabilities')
    parser.add_argument('--nfolds',
        type=int,
        default='none',
        help='Number of folds for cross validation')
    parser.add_argument('--confounders',
        type=str,
        help='Table with target dataset specific QC, e.g. SSMD')
    parser.add_argument('--feature-dir',
        type=str,
        help='Directory with the feature input files')
    parser.add_argument('--feature-info',
        type=str,
        help='Table feature datasets requiring dataset and filename columns')
    parser.add_argument('--model',
        type=str,
        help='The Name field for one of the models from the models config')
    parser.add_argument('--start-col',
        type=int,
        default=0,
        help='Start column index of Y matrix to determine range of genes to predict')
    parser.add_argument('--end-col',
        type=int,
        default=None,
        help='Number of folds for cross validation')
    parser.add_argument('--feat-suffix',
        type=str,
        default='features.csv',
        help='File for writing feature importance')
    parser.add_argument('--pred-suffix',
        type=str,
        default='predictions.csv',
        help='File for writing predictions')
    parser.add_argument('--output-dir',
        type=str,
        help='The output folder for storing the outputs')
    return parser



def soft_roc_auc(ytrue, ypred):
    if all(ytrue > .5) or not any(ytrue > .5):
        return .5
    return roc_auc_score((ytrue > .5).astype(np.int), ypred)

def pearson_score(ytrue, ypred):
    if all(ypred == np.mean(np.mean(ypred))) or all(ytrue == np.mean(np.mean(ytrue))):
        return 0
    return pearsonr(ytrue, ypred)[0]

def varfilter(X, threshold):
    return (np.std(np.array(X), axis=0) > threshold)

def pdnormalize(df):
    df[:] -= df.mean()
    df[:] /=df.std()

def align_df(df, suffix, cells):
    out = df.loc[cells]
    out = out[out.columns[(out != 0).any(axis=0)]]
    out.columns = [s+suffix for s in out.columns]
    return out

class QuantileKFold(KFold):
    """Quantile K-Folds cross-validator
    Sorts continuous values and then bins them into N/k number of bins,
    such that each bin contains k values.
    Values in each bin are then assigned a random fold 1:k
    """

    def __init__(self, n_splits=5, shuffle=False, random_state=None):
        super().__init__(n_splits=n_splits, shuffle=shuffle, random_state=random_state)


    def split(self, X, y, groups=None):
        """Generate indices to split data into training and test set.
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Always ignored, exists for compatibility.
        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.
            Stratification is done based on the y labels.
        groups : object
            Always ignored, exists for compatibility.
        Yields
        ------
        train : ndarray
            The training set indices for that split.
        test : ndarray
            The testing set indices for that split.
        Notes
        -----
        Randomized CV splitters may return different results for each call of
        split. You can make the results identical by setting ``random_state``
        to an integer.
        """
        #y = check_array(y, ensure_2d=False, dtype=None)
        df = pd.DataFrame(data=y)
        df.columns = ['value']
        df['position'] = range(0,df.shape[0])
        df=df.sort_values(by=['value'])

        fold_assignment = {}
        for i in range(0,math.ceil(df.shape[0]/self.n_splits)):
            list = [*range(0,self.n_splits)]
            random.shuffle(list)
            fold_assignment[i] = list
        fold_assignment = [*chain(*fold_assignment.values())]
        df['fold'] = [fold_assignment[i] for i in range(0,df.shape[0])]

        def split_pairs(df,nfolds):
            for i in range(nfolds):
                train_pos = df[df['fold'] != i]['position']
                test_pos = df[df['fold'] == i]['position']
                yield [train_pos.tolist(),test_pos.tolist()]

        split_generator = split_pairs(df=df,nfolds=self.n_splits)
        return split_generator


def single_fit(column, X, Y, model_types, splitter, scoring, nfeatures=10, rounding=False, return_models=False):
    y = Y[column]
    if rounding:
        if (y > .5).sum() == 0 or (y > .5).sum() == len(y):
            raise ValueError("Column %s has %i of %i possible true values\n%r" %(
                column, (y > .5).sum(), len(y), Y[column]
            )
                            )
    else:
        if y.std() == 0:
            raise ValueError("Column %s has 0 variance\n%r" %(
            column, y)
        )
    if y.isnull().any():
        raise ValueError('Column %s of y contains %i nulls' %(column, y.isnull().sum()))
    target_models = [val['ModelClass'](**val['kwargs']) for val in model_types]
    scores = []
    features = []
    prediction = []
    for model, x in zip(target_models, X):
        if x.isnull().any().any():
            raise ValueError('Feature set for model %r contains nulls. Axial sums of nulls:\n%r\n\n%r' %(
                model, x.isnull().sum(), x.isnull().sum(axis=1)))
        if rounding:
            splits = splitter.split(y > .5, y > .5)
        else:
            splits = splitter.split(y, y)
        score = []
        n = 0
        model_prediction = pd.Series(np.nan, index=Y.index, name=column)
        for train, test in splits:
            try:
                model.fit(x.iloc[train], y.iloc[train])
                ypred = model.predict(x.iloc[test])
            except Exception as e:
                print('error fitting model %r for column %s' %(model, column))
                print('train indices:\n %r\n' %train)
                print('test indices:\n %r\n' %test)
                print('train features: \n%r\n' %x.iloc[train])
                print('test features: \n%r\n' %x.iloc[test])
                print('train column: \n%r\n' %y.iloc[train])
                print('test column: \n%r\n' %y.iloc[test])
                raise e
            model_prediction.iloc[test] = ypred[:]
            score.append(scoring(y.iloc[test], ypred))
            n += 1
        scores.append(score)
        prediction.append(model_prediction)
        model.fit(x, y)
        try:
            features.append(model.get_feature_series(nfeatures))
        except AttributeError:
            features.append(pd.Series(dtype=np.float64))
        gc.collect()
    best_index = np.argmax(np.mean(scores, axis=1))
    if not return_models:
     target_models = [np.nan for i in range(len(scores))]
    return {'models': target_models, 'best': best_index, 'scores': scores, 'features': features, 'predictions': prediction}


##############################################################
#################### E N S E M B L E #########################
##############################################################

class EnsembleRegressor:
    def __init__(self, model_types, nfolds=3, scoring=soft_roc_auc, Splitter=StratifiedKFold, rounding=False):

        '''
        model_types: [{'Name': str, 'ModelClass': class,  'kwargs': dict}] Model classes will be initiated with the
            dict of keyword arguments
        '''
        self.model_types = model_types
        self.best_indices = {}
        self.trained_models = {}
        self.scores = {}
        self.important_features = {}
        self.nfolds = nfolds
        self.splitter = Splitter(n_splits=nfolds, shuffle=True)
        self.scoring = scoring
        self.columns = None
        self.rounding = rounding
        self.predictions = None

    def check_x(self, X):
        xerror = ValueError(
            'X must be a list or array with a feature set dataframe of matching indices for each model \
            present in the ensemble, passed\n%r'
            %X
        )
        if not len(X) == len(self.model_types):
            print('X not the same length as models\n')
            raise xerror
        for df in X[1:]:
            if not all(df.index == X[0].index):
                raise xerrorcheck_random_state

    def fit(self, X, Y, columns=None, report_freq=20):
        '''
        X: [{ModelClass: dataframe}
        Y: dataframe
        '''
        self.check_x(X)
        assert isinstance(Y, pd.DataFrame)
        if not all(Y.index == X[0].index):
            raise ValueError('Y must be a dataframe with index matching the indices in X')
        if columns is None:
            columns = Y.columns
        self.columns = Y.columns
        n = len(self.model_types)
        outputs = {'models': {}, 'best': {}, 'scores': {}, 'features': {}, 'predictions': {}}
        start_time = time()
        curr_time = start_time
        for i, col in enumerate(columns):
            ind = Y.index[Y[col].notnull()]
            output = single_fit(column=col, X=[x.loc[ind] for x in X], Y=Y.loc[ind], model_types=self.model_types, splitter=self.splitter,
                      scoring=self.scoring, rounding=self.rounding)
            for key in outputs.keys():
                    outputs[key][col] = output[key]
            t = time()
            if t - curr_time > report_freq:
                print(
                    '%f elapsed, %i%% complete, %f estimated remaining' %(
                        t - start_time, int(100*(i+1)/len(columns)), (t-start_time)*(len(columns)-i-1)*1./(i+1))
                    )
                curr_time = t
        self.trained_models.update(outputs['models'])
        self.best_indices.update(outputs['best'])
        self.scores.update(outputs['scores'])
        self.important_features.update(outputs['features'])
        predictions = [{col: val[j] for col, val in outputs['predictions'].items()} for j in range(n)]
        if self.predictions is None:
            self.predictions = [pd.DataFrame(v) for v in predictions]
        else:
            for i in range(len(self.model_types)):
                self.predictions[i] = self.predictions[i].join(outputs['predictions'][i])

    def predict(self, X):
        self.check_x(X)
        return pd.DataFrame(
            {column: self.trained_models[column][self.best_indices[column]].predict(X[self.best_indices[column]])
                            for column in self.columns},
                           index = X[0].index)

    def save_results(self, feat_outfile, pred_outfile):
        columns = ['gene', 'model']
        for i in range(self.nfolds):
            columns.append('score%i' %i)
        columns.append('best')
        for i in range(10):
            columns.extend(['feature%i' %i, 'feature%i_importance' %i])

        melted = pd.DataFrame(columns=columns)
        for gene in self.trained_models.keys():
            for i in range(len(self.model_types)):
                row = {
                    'gene': gene,
                    'model': self.model_types[i]['Name'],
                    'best': self.best_indices[gene] == i
                }
                for j in range(self.nfolds):
                    row['score%i' %j] = self.scores[gene][i][j]
                for j in range(10):
                    try:
                        row['feature%i' %j] = self.important_features[gene][i].index[j]
                        row['feature%i_importance' %j] = self.important_features[gene][i].iloc[j]
                    except IndexError:
                        row['feature%i' %j] = np.nan
                        row['feature%i_importance' %j] = np.nan
                melted = melted.append(row, ignore_index=True)
        melted.to_csv(feat_outfile, index=None)
        for model, pred in zip(self.model_types, self.predictions):
            pred.to_csv(pred_outfile,index_label="Row.name")


##############################################################
###################  O T H E R  ##############################
##############################################################

def gene_feature_filter(df, gene_name):
    genes = [x.split('_')[0].split(' ')[0] for x in df.columns]
    mask = np.array([s == gene_name for s in genes], dtype=bool)
    return mask


class SelfFeatureForest(RandomForestRegressor):
    '''
    Uses only features matching ("like") the target (+ `reserved_columns` which are always included).
    Uses the target's name with everything after the first space stripped off to match columns in the feature set.
    A custom target name can be passed with the call to "fit" to be used for finding related features.
    '''
    def __init__(self, reserved_columns=[], **kwargs):
        self.reserved_columns = reserved_columns
        RandomForestRegressor.__init__(self, **kwargs)
        self.feature_names = []

    def fit(self, X, y, name=None, **kwargs):
        if name is None:
            name = y.name
        self.name = name
        mask = X.columns.isin(self.reserved_columns)
        mask = mask | gene_feature_filter(X, self.name.split(' ')[0])
        self.feature_names = X.columns[mask].tolist()
        features = X.loc[:, self.feature_names]
        RandomForestRegressor.fit(self, features.values, y.loc[features.index].values, **kwargs)
        #RandomForestRegressor.fit(self, features.as_matrix(), y.loc[features.index].values, **kwargs)

    def predict(self, X, **kwargs):
        features = X.loc[:, self.feature_names]
        return RandomForestRegressor.predict(self, features.values, **kwargs)
        #return RandomForestRegressor.predict(self, features.as_matrix(), **kwargs)

    def get_feature_series(self, n_features=None):
        if n_features is None:
            n_features = len(self.feature_names)
        imp = pd.Series(self.feature_importances_, index=self.feature_names)
        return imp.sort_values(ascending=False)[:n_features]


class RelatedFeatureForest(SelfFeatureForest):
    '''
    Uses a two-column list of related features to select only features related to the target (+ `reserved_columns` which are always included).
    Uses the target's name with everything after the first space stripped off to match to related features.
    A custom target name can be passed with the call to "fit" to be used for finding related features.
    '''
    def __init__(self, reserved_columns=[], relations=pd.DataFrame(columns=['target', 'partner']), **kwargs):
        SelfFeatureForest.__init__(self, reserved_columns, **kwargs)
        self.relations = relations

    def fit(self, X, y, name=None, **kwargs):
        if name is None:
            name = y.name
        self.name = name
        self.related = self.relations.partner[self.relations.target == name.split(' ')[0]].tolist()
        self.related = list(set(self.related))
        if not self.name in self.related:
            self.related.append(self.name)
        mask = X.columns.isin(self.reserved_columns)
        for partner in self.related:
            mask = mask | gene_feature_filter(X, partner)
        self.feature_names = X.columns[mask]
        features = X[self.feature_names]
        print('Fitting related model with X shape', features.shape)
        RandomForestRegressor.fit(self, features, y.loc[features.index].values, **kwargs)


class KFilteredForest(RandomForestRegressor):
    '''
    Selects the top `k` features with highest correlation with the target and variance greater than `var_threshold`
    '''
    def __init__(self, k=1000, var_threshold=0, **kwargs):
        self.k = k
        RandomForestRegressor.__init__(self, **kwargs)
        self.filter = SelectKBest(score_func=f_regression, k=k)
        self.var_threshold = var_threshold

    def fit(self, X, y, **kwargs):
        self.mask1 = varfilter(X, self.var_threshold)
        x = X.loc[:, X.columns[self.mask1]]
        if x.shape[1] > self.k:
            self.filter.fit(x, np.array(y))
            self.mask2 = self.filter.get_support()
            x = x.loc[:, x.columns[self.mask2]]
        self.feature_names = x.columns.tolist()
        RandomForestRegressor.fit(self, x.values, y, **kwargs)
        #RandomForestRegressor.fit(self, x.as_matrix(), y, **kwargs)


    def predict(self, X, **kwargs):
        x = X.loc[:, self.feature_names]
        return RandomForestRegressor.predict(self, x.values, **kwargs)
        #return RandomForestRegressor.predict(self, x.as_matrix(), **kwargs)

    def get_feature_series(self, n_features):
        if n_features is None:
            n_features = len(self.feature_names)
        imp = pd.Series(self.feature_importances_, index=self.feature_names)
        return imp.sort_values(ascending=False)[:n_features]

class PandasForest(RandomForestRegressor):
    '''
    A simple wrapper for RandomForestRegressor that plays nice with dataframes and series instead of numpy arrays
    '''
    def fit(self, X, y, **kwargs):
        self.feature_names = X.columns.tolist()
        RandomForestRegressor.fit(self, X, y, **kwargs)

    def get_feature_series(self, n_features):
        if n_features is None:
            n_features = len(self.feature_names)
        imp = pd.Series(self.feature_importances_, index=self.feature_names)
        return imp.sort_values(ascending=False)[:n_features]


def assemble_feature_set(target_samples, matrices = {}, tables = {}, required = [], fill_na=True):
    """Column binds the feature matrices together into one table.

    Description.

    Args:
        matrices (dict): Key is dataset name and value is filepath to numeric matrix
        tables (dict): Key is dataset name and value is filepath to table
        target_samples: indeces of the Y target matrix to filter by
        required: list of datasets (by name) that are required for valid samples
            including the name of a table includes all the columns
        fill_na:

    Returns:
        dataframe: X

    Raises:
        AttributeError: The ``Raises`` section is a list of all exceptions
            that are relevant to the interface.
        ValueError: If `param2` is equal to `param1`.

    """

    features = {}
    valid_samples = set(target_samples)

    #Load matrices, filter indexes by target_samples, add matrix to features dict
    if len(matrices) > 0:
        for key, val in matrices.items():
            print(val)
            df = pd.read_csv(val).set_index('Row.name') #TODO use load function instead
            if key in required:
                valid_samples = valid_samples.intersection(set(df.index))
            df = df.loc[df.index.isin(valid_samples)]
            features.update({key: df})

    #Load tables, filter by target_samples, add each column as a feature in the dict
    if len(tables) > 0:
        for key, item in tables.items():
            df = pd.read_csv(item).set_index('Row.name') #TODO use load function instead
            if key in required:
                valid_samples = valid_samples.intersection(set(df.index))
            df = df.loc[df.index.isin(valid_samples)]
            #Determine which columns are numeric and add them as individual features as-is
            numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
            num_vars = df.select_dtypes(include=numerics)
            features.update({key: num_vars})
            #features.update({(s+"_"+key): num_vars.loc[:,num_vars.columns.isin([s])] for s in num_vars.columns})
            #Determine which columns are categorical and add each as a one-hot encoded matrix
            #TODO change the lineage info format when loading confounders from taiga so that the column name isn't
            #     already included in the value. This way the original column name can be appended to column names of
            #     one-hot matrix and it will protect against losing track of where features came from if 2 columns of
            #     table share the same value
            cat_vars = df.select_dtypes(exclude=numerics)
            for s in cat_vars.columns:
                onehot = pd.get_dummies(cat_vars[s])
                detail_name = s + "_" + key
                features.update({detail_name: onehot})

    #Get list of valid samples that are in any feature dataset
    if len(required) == 0:
        valid_samples = list(set.union(*tuple([set(features[key].index) for key in features.keys()])))
    else:
        valid_samples = list(valid_samples)
        #valid_samples = list(set.intersection(*tuple([set(features[key].index) for key in required_refined])))
    assert len(valid_samples) > 0, "No cell lines found that are present in all of the features %r" %required

    #Remove duplicate rows and filter for valid samples
    for key, df in list(features.items()):
        if sum(df.index.duplicated()):
            print('%s has %i duplicates' %(key, sum(df.index.duplicated())))
            features[key] = df[~df.index.duplicated()]
        features[key] = df.loc[df.index.isin(valid_samples)]

    #Zscore all non-binary columns
    for val in features.values():
        if (~val.fillna(0).isin([0, 1])).any().any(): #non-binary dataframe
            try:
                pdnormalize(val)
            except Exception as e:
                print(val)
                raise e
    if len(features) == 1:
        return list(features.values())[0]

    print('joining features')
    megafeature = pd.DataFrame(index=valid_samples)
    for key, val in features.items():
        try:
            out = val.reindex(valid_samples)    #adds rows of NaN if feature type is not required and is missing samples in valid_samples
            out = out[out.columns[(out != 0).any(axis=0)]]
            suffix = '_%s' %key
            out.columns = [s+suffix for s in out.columns]
            megafeature = megafeature.join(out)
        except KeyError:
            print(key)
            print(val[:5])
            print()
    megafeature = megafeature.astype(np.float64)
    megafeature.dropna(how='all', axis=0, inplace=True) #drops samples that have all missing values
    megafeature.dropna(how='all', axis=1, inplace=True) #drops variables that are all missing
    if fill_na:
        print('filling NaNs')
        megafeature.fillna(0, inplace=True) #THIS FILLS ALL NA with 0. Shouldn't it impute based on the average instead??????????
    return megafeature

def run_model(X, Y, model, nfolds, feat_output, pred_output, start_col=0, end_col=None, task='regress', relation_table=None):
    """Fit models for specified columns of Y using a selection of feature subsets from X.

    Function parameters should be documented in the ``Args`` section. The name
    of each parameter is required. The type and description of each parameter
    is optional, but should be included if not obvious.

    Args:
        param1 (int): The first parameter.
        param2 (:obj:`str`, optional): The second parameter. Defaults to None.
            Second line of description should be indented.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        bool: True if successful, False otherwise.

        The return type is optional and may be specified at the beginning of
        the ``Returns`` section followed by a colon.

        The ``Returns`` section may span multiple lines and paragraphs.
        Following lines should be indented to match the first line.

        The ``Returns`` section supports any reStructuredText formatting,
        including literal blocks::

            {
                'param1': param1,
                'param2': param2
            }

    Raises:
        AttributeError: The ``Raises`` section is a list of all exceptions
            that are relevant to the interface.
        ValueError: If `param2` is equal to `param1`.

    """

    if end_col is None:
        end_col = Y.shape[1]
    print('aligning features')
    shared_lines = list(set(X.index) & set(Y.index))
    assert len(shared_lines) > 0, "no shared lines found: \n\n features %r\n\n "
    Y = Y.loc[shared_lines]
    X = X.loc[shared_lines]
    print('Number of shared cell lines: ' + str(len(shared_lines)))


    if model['Exempt'] is not None:
        constant_features = [s for s in X.columns if any(s.endswith(end) for end in model['Exempt'])]
    else:
        constant_features = []

    new_model = {'Name': model['Name']}

    if (model['Relation'] == "All") & (X.shape[1] <= 1000):
        new_model['ModelClass'] = PandasForest
        new_model['kwargs'] = dict(max_depth=8, n_estimators=100, min_samples_leaf=5)
    if (model['Relation'] == "All") & (X.shape[1] > 1000):
        new_model['ModelClass'] = KFilteredForest
        new_model['kwargs'] = dict(max_depth=8, n_estimators=100, min_samples_leaf=5)
    elif model['Relation'] == "MatchTarget":
        new_model['ModelClass'] = SelfFeatureForest
        new_model['kwargs'] = dict(reserved_columns=constant_features, max_depth=8, n_estimators=100, min_samples_leaf=5)
    elif model['Relation'] == "MatchRelated":
        new_model['ModelClass'] = RelatedFeatureForest
        new_model['kwargs'] = dict(reserved_columns=constant_features, relations=relation_table[['target', 'partner']], max_depth=8, n_estimators=100, min_samples_leaf=5)

    models = [new_model]

    if len(X) != len(Y):
        raise RuntimeError('length of X and Y do not match (shapes %r and %r)' %(X.shape, Y.shape))
    Xtrain = [X]
    assert len(Xtrain) == len(models), "number of models %i does not match number of feature sets %i" %(len(models), len(Xtrain))
    for i, x in enumerate(Xtrain):
        assert x.shape[1] > 0, "feature set %i does not have any columns" %i
        assert all(x.index == Y.index), "feature set %i index does not match Y index\n\n%r" %(i, x.iloc[:5, :5])

    print('creating TDA ensemble')
    if task == 'classify':
        ensemble = EnsembleRegressor(model_types=models, nfolds=nfolds, rounding=True, Splitter=StratifiedKFold, scoring=soft_roc_auc)
    elif task == 'regress':
        ensemble = EnsembleRegressor(model_types=models, nfolds=nfolds, rounding=False, Splitter=QuantileKFold, scoring=r2_score)
    else:
        raise ValueError('task must be "classify" or "regress"')

    columns = Y.columns[start_col : end_col]
    ensemble.fit(Xtrain, Y, columns)
    print('Finished fitting')
    ensemble.save_results(feat_outfile=feat_output,pred_outfile=pred_output)

if __name__ == '__main__':
    args = getParser().parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    #Read the model definitions
    with open(args.model_config) as f:
        model_json = json.load(f)

    #Lookup the model input parameter by name in the model config
    model_names = list(map(lambda x: x['Name'], model_json))
    model_def = model_json[model_names.index(args.model)]

    #Task mode - classify or regress
    # task_mode = model_def["Task"]["mode"]
    task_mode = args.task_mode

    #Load the file-info table to get location of datasets
    file_info = pd.read_csv(args.feature_info,sep='\t')

    #Get file paths for required feature datasets
    file_info['filename'] = [os.path.join(args.feature_dir,s) for s in file_info['filename']]
    all_files = dict(zip(file_info['dataset'], file_info['filename']))
    if model_def['Features'] is not None:
        file_info = file_info.loc[file_info.dataset.isin(model_def['Features'])]
        feature_matrices = dict(zip(file_info['dataset'], file_info['filename']))
    else:
        feature_matrices = []

    feature_tables = {}
    if model_def['Confounders'] == 'True':
        feature_tables["Confounders"] = args.confounders

    #Load the target data and sample info file and create a cell line list
    Y = pd.read_csv(args.targets,index_col=0)
    if args.end_col is None:
        Y = Y.iloc[:,args.start_col:]
    else:
        Y = Y.iloc[:,args.start_col:args.end_col]
    cls = Y.index

    #Call the assemble features function with the feature matrix list, feature table list, CL list, and optional relation table
    X = assemble_feature_set(target_samples=cls, matrices=feature_matrices, tables=feature_tables, required=model_def['Required'])
    print('X dimensions: ' + str(X.shape[0]) + "," + str(X.shape[1]))

    #Intersect target cell lines with features
    Y = Y.loc[list(set(X.index) & set(Y.index))]
    print('Y dimensions after intersect X: ' + str(Y.shape[0]) + "," + str(Y.shape[1]))

    #Load a relation table if necessary
    if (model_def["Relation"] == "MatchTarget") | (model_def["Relation"] == "All"):
        related_table = None
    else:
        related_table = pd.read_csv(all_files[model_def["Relation"]])
        related_table.partner = [s.split(' ')[0] for s in related_table.partner]
        related_table.target = [s.split(' ')[0] for s in related_table.target]

    outfile_feat = args.output_dir + args.model + "_" + str(args.start_col) + "_" + str(args.end_col) + "_" + args.feat_suffix + '.csv'
    outfile_pred = args.output_dir + args.model + "_" + str(args.start_col) + "_" + str(args.end_col) + "_" + args.pred_suffix + '.csv'

    run_model(X, Y, model=model_def, nfolds=args.nfolds, feat_output=outfile_feat, pred_output=outfile_pred, task=task_mode, relation_table=related_table)
