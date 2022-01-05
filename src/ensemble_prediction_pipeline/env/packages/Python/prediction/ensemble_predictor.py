from __future__ import print_function
import numpy as np
import pandas as pd
from time import time
from sklearn.linear_model import LinearRegression, ElasticNet
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.metrics import roc_auc_score, make_scorer, r2_score
from sklearn.feature_selection import SelectKBest, VarianceThreshold, mutual_info_regression, f_regression
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from scipy.stats import pearsonr
import argparse
import os
#from taigapy import TaigaClient
import gc


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
            features.append(pd.Series(dtype=np.float))
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
    
    def save_results(self, name):
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
        melted.to_csv(name + '_summary.csv', index=None)
        for model, pred in zip(self.model_types, self.predictions):
            pred.to_csv('%s_%s_predictions.csv' %(name, model['Name']))


##############################################################
###################  O T H E R  ##############################
##############################################################
class LookUp(LinearRegression):
    def __init__(self):
        pass
    
    def process_x(self, X):
        make_dummy = False
        if len(X.shape) == 1:
            make_dummy = True
        elif X.shape[1] < 2:
            make_dummy = True
        if make_dummy:
            dummies = pd.get_dummies(X)
            return dummies
        else:
            return X
    
    def fit(self, X, y):
        '''
        X: 2D array or dataframe. Can either contain a single column of categorical labels, or be a multicolumn of one-hot encoded
        values
        y: 1D array or pandas Series
        '''
        x = self.process_x(X)
        self.feature_names = x.columns.tolist()
        self.means = np.array(y).dot(np.array(x)) * 1.0/np.sum(x, axis=0)
        self.means[pd.isnull(self.means)] = 0
        
    def predict(self, X):
        x = self.process_x(X)
        return x.dot(self.means)
    
    def get_feature_series(self, n_features=None):
        if n_features is None:
            n_features = len(self.feature_names)
        def importance(p):
            return np.sqrt((p-np.mean(p))**2)
        ents = importance(self.means)
        return pd.Series(
            ents/sum(ents), index=self.feature_names
                        ).sort_values(ascending=False)[:n_features]


def gene_feature_filter(df, gene_name):
    genes = [x.split('_')[0].split(' ')[0] for x in df.columns]
    mask = np.array([s == gene_name for s in genes], dtype=np.bool)
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


class KFilteredFeatureTypeForest(SelfFeatureForest, KFilteredForest):
    '''
    Selects the top `k` highest correlation features out of a subset of features 
    based on their suffix (text occuring after the final underscore in the feature 
    column name). 
    '''
    def __init__(self, reserved_columns=[], suffixes=[], k=1000, var_threshold=0, **kwargs):
        SelfFeatureForest.__init__(self, reserved_columns, **kwargs)
        self.k = k
        self.filter = SelectKBest(score_func=f_regression, k=k)
        self.var_threshold = var_threshold
        self.suffixes = suffixes

    def fit(self, X, y, **kwargs):
        mask = X.columns.to_series().apply(lambda x: str(x).split('_')[-1] in self.suffixes)
        if mask.sum() < 1:
            unique_suffixes = sorted(set([str(s).split('_')[-1] for s in X.columns]))
            raise ValueError('None of the given suffixes %r were found in any of the feature suffixes %r' %(
                self.suffixes, unique_suffixes))
        KFilteredForest.fit(self, X[X.columns[mask]], y)


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

