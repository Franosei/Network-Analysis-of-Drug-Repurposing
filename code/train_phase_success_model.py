import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score
)

# Load dataset
df = pd.read_csv("graph/training_dataset.csv")

# Features and target
X = df[["DegreeCentrality", "EigenvectorCentrality", "BetweennessCentrality","ClusteringCoefficient", "RandomWalkCentrality"]]
y = df["PhaseSuccess"]

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42, stratify=y
)

# Scale features for models that benefit from normalization
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Define models
models = {
    "Logistic Regression (Ridge)": LogisticRegressionCV(
        cv=5, penalty="l2", solver="lbfgs", max_iter=1000
    ),
    "Logistic Regression (Lasso)": LogisticRegressionCV(
        cv=5, penalty="l1", solver="liblinear", max_iter=1000
    ),
    "Random Forest": RandomForestClassifier(n_estimators=100, random_state=42),
    "XGBoost": XGBClassifier(use_label_encoder=False, eval_metric="logloss", random_state=42)
}

# Train and evaluate
for name, model in models.items():
    print(f"\n==== {name} ====")
    if "Logistic Regression" in name:
        model.fit(X_train_scaled, y_train)
        y_pred = model.predict(X_test_scaled)
    else:
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

    # Classification report
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred, digits=4))

    # Confusion matrix
    print("Confusion Matrix:")
    print(confusion_matrix(y_test, y_pred))

    # Metrics
    acc = accuracy_score(y_test, y_pred)
    prec = precision_score(y_test, y_pred)
    rec = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)

    print(f"Accuracy:  {acc:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall:    {rec:.4f}")
    print(f"F1 Score:  {f1:.4f}")

    # Coefficients or Feature Importances
    if "Logistic Regression" in name:
        print("\nModel Coefficients (α):")
        for feature, coef in zip(X.columns, model.coef_[0]):
            print(f"  {feature}: {coef:.4f}")
    elif name == "Random Forest":
        print("\nFeature Importances:")
        for feature, imp in zip(X.columns, model.feature_importances_):
            print(f"  {feature}: {imp:.4f}")
    elif name == "XGBoost":
        print("\nFeature Importances:")
        for feature, imp in zip(X.columns, model.feature_importances_):
            print(f"  {feature}: {imp:.4f}")
