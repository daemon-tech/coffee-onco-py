# 001 – TCGA BRCA Subtypes Classification

**Goal**: Classify breast cancer subtypes (LumA, LumB, Basal, HER2) from TCGA gene expression data using ML.

**Datasets**: TCGA-BRCA (RNA-seq + clinical), downloaded 2026-01 via GDC API.

### Load Data from GDC API

**Option 1: Using the script**
```bash
cd projects/001_tcga_brca_subtypes
python load_data.py
```

**Option 2: Using Jupyter Notebook**
```bash
cd projects/001_tcga_brca_subtypes
jupyter notebook notebooks/01_load_tcga_data.ipynb
```

PAM50 subtypes are not directly available via the GDC API.