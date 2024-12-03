# RNA Editing Multiomic Machine Learning Models for Breast Cancer

This repository contains the code and resources necessary to reproduce analyses conducted as part of a research project in which we integrate multi-omics data with a particular focus on A-to-I(G) RNA editing to predict drug response in BC using machine learning (ML) models. Additionally, we designed a risk score for non-response based on ML-derived features using clinical trial data. While RNA editing events associated with drug response have been explored in lung cancer, gastric cancer, and AML, in BC, they have only been linked to survival outcomes. To address this gap, we also developed a simple probability-based score to assess therapy response. 

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
