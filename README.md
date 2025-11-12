# BOSS: Brain Orchestra Synchronization Study

BOSS (Brain Orchestra Synchronization Study) is a modular **MATLAB toolbox** for analyzing **feature loadings in synchronization activity** between brain regions using fMRI data.

It supports:
- **Canonical correlation analysis (CCA)** in a searchlight framework  
- **Mean-to-mean correlation** between regions  
- Flexible **ROI / atlas** handling and whole-brain searchlights  
- **Feature-based** connectivity (acoustic & emotion features, residual models)  

> Repository: https://github.com/mlionello/BOSS

---

## Key Features

- **Two-stage analysis pipeline**
  - **Group-level preparation**
    - Build common functional masks across subjects  
    - Compute group-average functional volumes  
    - Initialize models and sphere centers for searchlights  
  - **Subject-level analysis**
    - Apply group-fitted models to each subject  
    - Compute full and residual connectivity maps  
    - Aggregate subject-level results for group inference  

- **Flexible connectivity metrics**
  - Canonical correlation (CCA) between source and target regions  
  - CCA with dimensionality reduction (PCA on source or target)  
  - Simple correlation between mean signals of source and target  

- **Multiple analysis strategies**
  - Source region vs **all atlas regions**  
  - Source vs **single target region**  
  - Whole-brain **searchlight** connectivity  
  - Searchlight constrained within a **single ROI**  

- **Statistically robust inference**
  - Permutation-based **max-statistic** correction at the voxel/sphere level  
  - Paired t-tests between **full** and **residual** models  
  - Family-wise error–controlled significance maps  

- **Rich feature integration**
  - **Acoustic features** (MIR Toolbox, Essentia, VGGish)  
  - **Emotion features** (continuous slider ratings, deep emotion models)  
  - Joint PCA to build compact feature spaces  

---

## Requirements

- MATLAB (object-oriented features and standard toolboxes)
- fMRI data with derivatives organized in a directory referenced in:
  - `boss/constants.m`
- Optional/external (for full feature pipelines):
  - **MIR Toolbox** (low/mid-level acoustic features)
  - **Essentia** (spectral/timbral/musicological descriptors)
  - **VGGish** (deep audio embeddings)
  - AFNI (for clusterization and meta-analytic decoding, optional)
  - NiMARE / Neurosynth data (for meta-analytic decoding, optional)

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/mlionello/BOSS.git
