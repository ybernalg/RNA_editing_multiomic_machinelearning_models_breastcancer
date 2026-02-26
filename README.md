# RNA Editing Multi-omic Machine Learning Models for Breast Cancer

This repository contains code and resources to reproduce analyses from a research project focused on **A>I(G) RNA editing** and its contribution to predicting response to genotoxic therapies in breast cancer using **multi-omics** and **machine learning**.

---

## Project Description

In this study, we integrate multi-omics data—placing special emphasis on **A>I(G) RNA editing**—to predict treatment response in breast cancer using supervised machine learning models. The workflow combines RNA editing features with gene expression, DNA-derived features, and clinical covariates to evaluate whether RNA editing improves predictive performance.

---

## Key Components of the Project

### 1. Preprocessing
Data preparation, characterization, and harmonization across:
- Gene expression (EXP)
- RNA editing sites (ED)
- Tumoral and germline DNA features (DNA)
- Clinical variables (one-hot encoded)

Key preprocessing steps include:
- Mapping sequencing run IDs to subject/sample identifiers
- Aligning samples shared across omics layers and clinical data
- Removing zero-variance features (training set only)
- Train-only differential filtering for RNA editing and gene expression features
- Prefix-based feature tracking (ED__/EXP__/DNA)

All feature selection steps are performed **exclusively within the training set** to prevent data leakage.

---

### 2. Machine Learning Models

Classification models implemented within a unified training/testing framework:

- **Generalized Linear Models (GLM; elastic net via glmnet)**
- **Random Forest**
- **Support Vector Machines (linear SVM)**

Feature selection includes:
- **LASSO (glmnet)** for sparse model selection

If too few omics features remain after filtering, the pipeline automatically falls back to clinical-only models.

---

### 3. Model Evaluation

Models are evaluated on held-out test sets across multiple random seeds.

Primary performance metrics:

- **F1-score** (positive class = non-responder, `"no"`)
- **Matthews Correlation Coefficient (MCC)**
- **Precision–Recall Area Under the Curve (PR-AUC)**

Performance is summarized across seeds using:
- Mean
- Standard deviation
- Percentile-based 95% confidence intervals

---

---
