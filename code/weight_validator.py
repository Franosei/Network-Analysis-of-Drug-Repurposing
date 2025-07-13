import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.utils import resample
from sklearn.metrics import accuracy_score, roc_auc_score


class LogisticWeightValidator:
    """
    Validate the robustness of logistic regression feature weights using:
    1. Bootstrap resampling
    2. Gaussian noise injection
    3. Feature ablation

    This is designed to ensure that the learned weights are stable and interpretable,
    especially in highly separable datasets.
    """

    def __init__(self, X: pd.DataFrame, y: pd.Series, feature_names: list):
        """
        Parameters:
        - X: Feature matrix (DataFrame)
        - y: Binary labels (Series or ndarray)
        - feature_names: List of feature names in X
        """
        self.X = X.copy()
        self.y = y.copy()
        self.feature_names = feature_names
        self.model = LogisticRegression(max_iter=1000, class_weight='balanced')

    def bootstrap_weights(self, n_bootstrap: int = 100, random_state: int = 42) -> pd.DataFrame:
        """
        Fit logistic regression on multiple bootstrap resamples of the training data.
        Returns a DataFrame of weight estimates for each feature across bootstrap runs.
        """
        np.random.seed(random_state)
        weights = []

        for _ in range(n_bootstrap):
            X_sample, y_sample = resample(self.X, self.y)
            self.model.fit(X_sample, y_sample)
            weights.append(self.model.coef_[0])

        weights_df = pd.DataFrame(weights, columns=self.feature_names)
        return weights_df

    def noise_injection_test(self, noise_std: float = 0.01, n_runs: int = 50, random_state: int = 42) -> pd.DataFrame:
        """
        Add small Gaussian noise to the features and retrain to assess weight stability.
        Returns a DataFrame of weights over multiple noisy runs.
        """
        np.random.seed(random_state)
        weights = []

        for _ in range(n_runs):
            noise = np.random.normal(0, noise_std, size=self.X.shape)
            X_noisy = self.X.values + noise
            self.model.fit(X_noisy, self.y)
            weights.append(self.model.coef_[0])

        weights_df = pd.DataFrame(weights, columns=self.feature_names)
        return weights_df

    def ablation_test(self) -> pd.DataFrame:
        """
        Perform feature ablation by removing one feature at a time and retraining.
        Records test performance and weight changes.
        Returns a DataFrame summarizing performance drop and weight impact.
        """
        baseline_model = LogisticRegression(max_iter=1000, class_weight='balanced')
        baseline_model.fit(self.X, self.y)
        baseline_weights = baseline_model.coef_[0]
        baseline_acc = accuracy_score(self.y, baseline_model.predict(self.X))
        baseline_auc = roc_auc_score(self.y, baseline_model.predict_proba(self.X)[:, 1])

        results = []

        for idx, feature in enumerate(self.feature_names):
            X_reduced = self.X.drop(columns=[feature])
            model = LogisticRegression(max_iter=1000, class_weight='balanced')
            model.fit(X_reduced, self.y)
            acc = accuracy_score(self.y, model.predict(X_reduced))
            auc = roc_auc_score(self.y, model.predict_proba(X_reduced)[:, 1])
            delta_acc = baseline_acc - acc
            delta_auc = baseline_auc - auc
            results.append({
                "RemovedFeature": feature,
                "DeltaAccuracy": delta_acc,
                "DeltaAUC": delta_auc,
                "WeightBeforeRemoval": baseline_weights[idx]
            })

        return pd.DataFrame(results)

    def summary(self):
        """
        Print summary of bootstrap and noise test results for interpretation.
        """
        print("\n--- Bootstrap Weight Stability ---")
        boot_df = self.bootstrap_weights()
        print(boot_df.describe())

        print("\n--- Noise Injection Weight Stability ---")
        noise_df = self.noise_injection_test()
        print(noise_df.describe())

        print("\n--- Feature Ablation Test ---")
        ablation_df = self.ablation_test()
        print(ablation_df.to_string(index=False))


# === MAIN ===

if __name__ == "__main__":
    # Load separate known and unknown datasets
    known_path = "graph/graph_features_known.csv"
    unknown_path = "graph/graph_features_unknown.csv"

    known = pd.read_csv(known_path)
    unknown = pd.read_csv(unknown_path)

    # Assign labels
    known["Label"] = 1
    unknown["Label"] = 0

    # Combine datasets
    data = pd.concat([known, unknown], ignore_index=True)

    # Define features
    feature_names = [
        "GraphDistanceToIndication",
        "RandomWalkScore",
        "StructuralLikelihood",
        "PreferentialAttachment",
        "KatzSimilarity"
    ]

    X = data[feature_names]
    y = data["Label"]

    validator = LogisticWeightValidator(X, y, feature_names)
    validator.summary()
