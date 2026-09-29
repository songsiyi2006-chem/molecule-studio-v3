"""Small, serially fitted validation-selected ensembles with inspectable members."""
import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin, TransformerMixin, clone
from sklearn.utils.validation import check_is_fitted


class FeatureSubset(TransformerMixin, BaseEstimator):
    """Select a fixed representation without mutating shared float32 input."""
    def __init__(self, indices):
        self.indices = indices

    def fit(self, X, y=None):
        self.n_features_in_ = np.asarray(X).shape[1]
        return self

    def transform(self, X):
        check_is_fitted(self, "n_features_in_")
        values = np.asarray(X)
        if values.ndim != 2 or values.shape[1] != self.n_features_in_:
            raise ValueError("Feature schema mismatch")
        return np.array(values[:, self.indices], dtype=np.float32, order="C", copy=True)


class ValidationEnsemble(RegressorMixin, BaseEstimator):
    """One or two predefined estimators, fixed nonnegative validation weights.

Weight selection occurs in the training script on validation predictions only.
Fit never adjusts weights. Member disagreement is not a confidence interval.
"""
    def __init__(self, members, weights):
        self.members = members
        self.weights = weights

    def fit(self, X, y):
        weights = np.asarray(self.weights, dtype=float)
        if not 1 <= len(self.members) <= 2 or len(weights) != len(self.members):
            raise ValueError("Expected one or two weighted members")
        if not np.isfinite(weights).all() or (weights < 0).any() or not np.isclose(weights.sum(), 1):
            raise ValueError("Weights must be finite, nonnegative and sum to one")
        self.n_features_in_ = np.asarray(X).shape[1]
        self.member_names_ = [name for name, _ in self.members]
        self.member_weights_ = weights
        self.estimators_ = [(name, clone(estimator).fit(X, y)) for name, estimator in self.members]
        return self

    def predict_members(self, X):
        check_is_fitted(self, "estimators_")
        if np.asarray(X).ndim != 2 or np.asarray(X).shape[1] != self.n_features_in_:
            raise ValueError("Feature schema mismatch")
        return np.column_stack([estimator.predict(X) for _, estimator in self.estimators_])

    def predict(self, X):
        return self.predict_members(X) @ self.member_weights_
