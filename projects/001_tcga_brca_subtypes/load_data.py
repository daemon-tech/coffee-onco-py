# Script to load TCGA-BRCA data from GDC API.

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from data_loader import TCGADataLoader

def main():
    # Initialize loader
    loader = TCGADataLoader(data_dir="data")
    
    # Load all data (with auto-download of PAM50 if not found)
    manifest, clinical, pam50 = loader.load_all(auto_download_pam50=True)
    
    # Print summary
    print("\n" + "=" * 60)
    print("Data Summary")
    print("=" * 60)
    print(f"\nFile Manifest:")
    print(f"  Total files: {len(manifest)}")
    print(f"  Columns: {list(manifest.columns)}")
    
    print(f"\nClinical Data:")
    print(f"  Total cases: {len(clinical)}")
    print(f"  Columns: {list(clinical.columns)}")
    
    if not pam50.empty:
        print(f"\nPAM50 Subtypes:")
        print(f"  Total annotations: {len(pam50)}")
        print(f"  Columns: {list(pam50.columns)}")
        if 'pam50_subtype' in pam50.columns:
            print(f"  Subtype distribution:")
            print(pam50['pam50_subtype'].value_counts().to_string())
    else:
        print(f"\nPAM50 Subtypes: Not found")
        print(f"\n  Try manually downloading:")
        print(f"    file_path = loader.download_pam50_from_url()")
        print(f"    pam50 = loader.load_pam50_from_file(str(file_path))")
    
    print("\n" + "=" * 60)
    print("Next steps:")
    print("1. Review the data in data/ directory")
    print("2. Download expression files using loader.download_file()")
    print("3. Parse expression data and merge with clinical annotations")
    print("4. Use PAM50 subtypes for classification modeling")
    print("=" * 60)


if __name__ == "__main__":
    main()
