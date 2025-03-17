# Functions for the analyses related to TNBC and circadian rythmicity
import pandas as pd
import numpy as np

from sklearn.linear_model import LinearRegression, ElasticNet, LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import LeaveOneOut, GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import pearsonr


def get_correlations(df, y):
    """
    Calculates the correlations between each column of a DataFrame and a Series.
    To be used when calculating the correlations of each feature to a target.
    :param df: Pandas DataFrame with N columns and P rows
    :param y: Pandas Series with P entries
    :return: a Pandas Series containing the correlation coefficient (Pearson) for each column
    """
    if df is not None and y is not None:
        correls = pd.Series(name='correlation')
        pvals = pd.Series(name='p.values')
        for gene in df.columns:
            data = pd.concat([df[gene], y], axis=1)
            correls[gene] = (abs(data.corr().iloc[0, -1]))
            pvals[gene] = (pearsonr(data[gene], y)[1])
        return correls, pvals
    else:
        raise ValueError('Either df or y is not set')


def get_model(X, y, X_test, y_test, model_type):
    MODEL_TYPES = {'regression_forest': RandomForestRegressor(max_depth=10,
                                                              n_estimators=100,
                                                              criterion='squared_error',
                                                              max_features=0.5,
                                                              random_state=42, ),
                   'linear_regression': LinearRegression(),
                   'elastic_net': ElasticNet(random_state=42),
                   'classification_forest': RandomForestClassifier(max_depth=10, max_features=0.5, random_state=42),
                   'classification_tree': DecisionTreeClassifier(max_depth=10, random_state=42),
                   'nearest-neighbor': KNeighborsClassifier(n_neighbors=3),
                   'logistic': LogisticRegression(fit_intercept=False, random_state=42, max_iter=1000),
                   }
    if X is not None and y is not None:
        m = None
        try:
            m = MODEL_TYPES[model_type]
        except ValueError:
            print(f'model type {model_type} not in {MODEL_TYPES.keys()}')
        try:
            m.fit(X, y)
        except:
            print(f'not possible to fit the model to the data')
        try:
            prediction = m.predict(X_test.to_frame().transpose()).squeeze()
        except:
            try:
                pass  # prediction = m.predict(X_test.to_numpy().reshape(1, -1))
            except:
                print(f'could not form a prediction')
        try:
            val = (mean_squared_error(y_test, prediction)) ** 0.5
            print(f'RMSE: {val}')
        except:
            if y_test == prediction:
                val = 1
            else:
                val = 0
        return m, prediction, val
    else:
        raise ValueError('Either df or y is not set')


def grid_search_loocv(X, y, model, param_grid):
    loocv = LeaveOneOut()
    gs = GridSearchCV(model, param_grid, scoring='neg_mean_squared_error', cv=loocv, verbose=1)
    gs.fit(X, y)
    RMSE = pd.DataFrame(gs.cv_results_).filter(regex='split(\d|10)_')  # currently only works for 10 splits max
    RMSE['mean_test_RMSE'] = RMSE.abs().mean(axis=1) ** 0.5  # mean across rows
    return pd.concat([pd.DataFrame(gs.cv_results_), RMSE['mean_test_RMSE']], axis=1)
