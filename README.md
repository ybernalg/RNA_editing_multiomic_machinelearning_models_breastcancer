# RNA Editing Multiomic Machine Learning Models for Breast Cancer

This repository contains the code and resources necessary to reproduce analyses conducted as part of a research project focusing on **A>I(G) RNA editing** and its impact on the response to genotoxic drugs in breast cancer. The study leverages multi-omics approaches and machine learning models to uncover key biological insights.


---

## Project Description

In this study, we integrate multi-omics data with a particular focus on **A>I(G) RNA editing to predict drug response in BC using machine learning (ML) models**. Additionally, we designed a risk score for non-response based on ML-derived features using clinical trial data. While RNA editing events associated with drug response have been explored in lung cancer, gastric cancer, and AML, in BC, they have only been linked to survival outcomes. To address this gap, we also developed a simple probability-based score to assess therapy response. 

---

## Key Components of the Project

### 1. **Preprocessing**
   - Data preparation and visualization for:
     - Gene expression
     - RNA edited sites
     - Tumoral and germline DNA

### 2. **Machine Learning Models**
   - Implementation of various classification models, including:
     - **Generalized Linear Models (GLM)**
     - **Random Forest**
     - **Support Vector Machines (SVM)**
   - Feature selection strategies using:
     - **LASSO**
     - **Principal Component Analysis (PCA)**
   - Performance evaluation through:
     - **F1-score**
     - **AUC ROC curve**

### 3. **Risk Scoring**
   - Development of a **risk score** for therapy non-response based on features derived from the final machine learning model.

---
