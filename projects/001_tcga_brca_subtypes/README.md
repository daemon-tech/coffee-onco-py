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

**Option 3: Using Python interactively**
```python
from src.data_loader import TCGADataLoader

loader = TCGADataLoader(data_dir="data")
manifest, clinical, pam50 = loader.load_all()
```

## Data Loading

The `TCGADataLoader` class provides methods to:
- Query GDC API for TCGA-BRCA files
- Download file manifests for RNA-seq expression data
- Fetch clinical annotations
- Attempt to retrieve PAM50 subtype labels
- Download individual files by file ID

Note: PAM50 subtypes may need to be obtained from supplementary files or TCGA publications if not available via the annotations API.

## PAM50 Subtypes

PAM50 subtypes are not directly available via the GDC API. The `get_pam50_subtypes()` method tries multiple approaches:

1. **Annotations endpoint** - Checks for PAM50 annotations
2. **Supplementary files** - Searches for files with "PAM50" or "subtype" in the name
3. **Clinical data** - Checks if subtype information is in clinical fields

If PAM50 subtypes are not found, you can:

### Option 1: Load from downloaded file
```python
loader = TCGADataLoader(data_dir="data")
pam50 = loader.load_pam50_from_file("path/to/BRCA.547.PAM50.SigClust.Subtypes.txt")
```

### Option 2: Download from TCGA publications
- **Nature 2012 paper**: https://gdc.cancer.gov/about-data/publications/brca_2012
- File: `BRCA.547.PAM50.SigClust.Subtypes.txt`
- Contains PAM50 assignments for 547 samples

### Option 3: Use TCGAbiolinks (R package)
```r
library(TCGAbiolinks)
subtypes <- TCGAquery_subtype(tumor = "BRCA")
# Export to CSV and load in Python
```

### Option 4: Infer from expression data
Use PAM50 gene signature (50 genes) to classify samples based on expression patterns. This requires:
- Normalized expression data
- PAM50 gene list
- Classification algorithm (centroid-based or similar)