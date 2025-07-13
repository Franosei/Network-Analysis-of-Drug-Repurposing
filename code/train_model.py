import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report


class GraphFeatureWeightLearner:
    """
    Learn interpretable feature weights from known and unknown drug-disease graph features
    using various models (logistic, ridge, l1-penalized logistic, random forest, XGBoost).
    Evaluates models and reports weights or importances for Bayesian likelihood construction.
    Handles class imbalance using automatic class weighting. Avoids data leakage by isolating preprocessing to training data.
    """

    def __init__(self, known_path, unknown_path):
        """
        Parameters:
        - known_path: path to known (label=1) drug-disease features CSV
        - unknown_path: path to unknown (label=0) drug-disease features CSV
        """
        self.known_path = known_path
        self.unknown_path = unknown_path
        self.feature_names = [
            "GraphDistanceToIndication",
            "RandomWalkScore",
            "StructuralLikelihood",
            "PreferentialAttachment",
            "KatzSimilarity"
        ]
        self.models = {}

    def load_and_prepare_data(self):
        """
        Load, combine, and preprocess the known and unknown datasets.
        Standardizes features using statistics computed only on the training set to avoid data leakage.
        """
        known = pd.read_csv(self.known_path)
        unknown = pd.read_csv(self.unknown_path)

        known["Label"] = 1
        unknown["Label"] = 0

        combined = pd.concat([known, unknown], ignore_index=True)
        X = combined[self.feature_names]
        y = combined["Label"]

        # Split before scaling to avoid leakage
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            X, y, test_size=0.3, random_state=42, stratify=y
        )

        # Standardize only using training set
        self.scaler = StandardScaler()
        self.X_train = self.scaler.fit_transform(self.X_train)
        self.X_test = self.scaler.transform(self.X_test)

    def train_models(self):
        """
        Train and evaluate five classifiers.
        Applies class weighting only to training. Test set remains untouched.
        """
        self.models = {
            "LogisticRegression": LogisticRegression(max_iter=1000, class_weight="balanced"),
            "RidgeClassifier": RidgeClassifier(class_weight="balanced"),
            "LogisticL1": LogisticRegression(penalty="l1", solver="liblinear", max_iter=1000, class_weight="balanced"),
            "RandomForest": RandomForestClassifier(n_estimators=100, random_state=42, class_weight="balanced"),
            "XGBoost": XGBClassifier(use_label_encoder=False, eval_metric="logloss", random_state=42, scale_pos_weight=self._get_xgb_pos_weight())
        }

        for name, model in self.models.items():
            model.fit(self.X_train, self.y_train)
            y_pred = model.predict(self.X_test)
            y_prob = model.predict_proba(self.X_test)[:, 1] if hasattr(model, "predict_proba") else y_pred

            print(f"\n=== {name} ===")
            print("Accuracy:", accuracy_score(self.y_test, y_pred))
            print("ROC AUC:", roc_auc_score(self.y_test, y_prob))
            print(classification_report(self.y_test, y_pred))

    def _get_xgb_pos_weight(self):
        """
        Compute scale_pos_weight for XGBoost based on class distribution in training set.
        This balances the loss for imbalanced classes.
        """
        n_pos = sum(self.y_train == 1)
        n_neg = sum(self.y_train == 0)
        return n_neg / n_pos if n_pos > 0 else 1.0

    def extract_weights(self):
        """
        Extract feature weights or importances from trained models.
        Linear models expose coefficients, tree models expose feature importances.
        """
        print("\n--- Feature Weights / Importances ---")
        for name, model in self.models.items():
            print(f"\n{name}")
            if name in ["LogisticRegression", "LogisticL1"]:
                weights = model.coef_[0]
                for f, w in zip(self.feature_names, weights):
                    print(f"{f:30s}: {round(w, 4)}")

            elif name == "RidgeClassifier":
                weights = model.coef_
                for f, w in zip(self.feature_names, weights):
                    print(f"{f:30s}: {round(w, 4)}")

            elif name in ["RandomForest", "XGBoost"]:
                importances = model.feature_importances_
                for f, w in zip(self.feature_names, importances):
                    print(f"{f:30s}: {round(w, 4)}")

    def run(self):
        """
        Full pipeline: load data, train models, evaluate, and extract interpretable results.
        """
        self.load_and_prepare_data()
        self.train_models()
        self.extract_weights()


# ==== Run ====
if __name__ == "__main__":
    learner = GraphFeatureWeightLearner(
        known_path="graph/graph_features_known.csv",
        unknown_path="graph/graph_features_unknown.csv"
    )
    learner.run()
