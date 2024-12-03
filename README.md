# RNA Editing Multiomic Machine Learning Models for Breast Cancer

This repository contains the code and resources necessary to reproduce analyses conducted as part of a research project focusing on **A>I(G) RNA editing** and its impact on the response to genotoxic drugs in breast cancer. The study leverages multi-omics approaches and machine learning models to uncover key biological insights.

---

## Project Description

The primary objective of this project is to investigate the relationship between **RNA editing** and **therapy response** in breast cancer patients. By integrating **transcriptomic**, **genomic**, and **clinical data**, machine learning models are employed to identify RNA editing patterns potentially associated with treatment resistance or sensitivity.

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
