# TCGA-BRCA Data Loader from GDC (Genomic Data Commons) API. This module handles downloading and loading TCGA-BRCA RNA-seq and clinical data from the GDC Data Portal API.

import os
import json
import requests
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Tuple
import time


class TCGADataLoader:
    # Load TCGA-BRCA data from GDC API.
    
    GDC_API_BASE = "https://api.gdc.cancer.gov"
    PROJECT_ID = "TCGA-BRCA"
    
    def __init__(self, data_dir: str = "data"):
        # Initialize the data loader. Parameters: data_dir (str) - Directory to save downloaded data files
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        
    def _query_gdc(self, endpoint: str, params: Dict, max_retries: int = 3) -> Dict:
        # Query GDC API with retry logic. Parameters: endpoint (str) - API endpoint (e.g., '/files', '/cases'), params (dict) - Query parameters, max_retries (int) - Maximum number of retry attempts. Returns: dict - JSON response from API
        url = f"{self.GDC_API_BASE}{endpoint}"
        
        for attempt in range(max_retries):
            try:
                response = self.session.get(url, params=params, timeout=60)
                response.raise_for_status()
                return response.json()
            except requests.exceptions.RequestException as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"Request failed, retrying in {wait_time}s... ({e})")
                    time.sleep(wait_time)
                else:
                    raise
    
    def get_file_manifest(self, data_type: str = "Gene Expression Quantification") -> pd.DataFrame:
        # Get file manifest for TCGA-BRCA RNA-seq expression data. Parameters: data_type (str) - Type of data to query (default: "Gene Expression Quantification"). Returns: pd.DataFrame - DataFrame with file IDs, names, and metadata
        print(f"Querying GDC API for {data_type} files...")
        
        # Query for files
        files_endpoint = "/files"
        params = {
            "filters": json.dumps({
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "cases.project.project_id",
                            "value": [self.PROJECT_ID]
                        }
                    },
                    {
                        "op": "in",
                        "content": {
                            "field": "files.data_type",
                            "value": [data_type]
                        }
                    },
                    {
                        "op": "in",
                        "content": {
                            "field": "files.experimental_strategy",
                            "value": ["RNA-Seq"]
                        }
                    }
                ]
            }),
            "fields": "file_id,file_name,file_size,cases.case_id,cases.samples.sample_id",
            "size": "10000",  # Adjust if needed
            "format": "JSON"
        }
        
        response = self._query_gdc(files_endpoint, params)
        
        if "data" not in response or "hits" not in response["data"]:
            raise ValueError("No files found in GDC response")
        
        files = response["data"]["hits"]
        print(f"Found {len(files)} files")
        
        # Flatten the nested structure
        file_list = []
        for file_info in files:
            file_list.append({
                "file_id": file_info.get("id"),
                "file_name": file_info.get("file_name"),
                "file_size": file_info.get("file_size"),
                "case_id": file_info.get("cases", [{}])[0].get("case_id") if file_info.get("cases") else None,
                "sample_id": file_info.get("cases", [{}])[0].get("samples", [{}])[0].get("sample_id") 
                            if file_info.get("cases") and file_info.get("cases")[0].get("samples") else None
            })
        
        manifest_df = pd.DataFrame(file_list)
        manifest_path = self.data_dir / "file_manifest.tsv"
        manifest_df.to_csv(manifest_path, sep="\t", index=False)
        print(f"Saved manifest to {manifest_path}")
        
        return manifest_df
    
    def get_clinical_data(self) -> pd.DataFrame:
        # Get clinical data for TCGA-BRCA cases. Returns: pd.DataFrame - DataFrame with clinical annotations including PAM50 subtypes
        print("Querying GDC API for clinical data...")
        
        # Query for cases with clinical data
        cases_endpoint = "/cases"
        params = {
            "filters": json.dumps({
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [self.PROJECT_ID]
                }
            }),
            "fields": "case_id,demographic,diagnoses,exposures",
            "size": "10000",
            "format": "JSON"
        }
        
        response = self._query_gdc(cases_endpoint, params)
        
        if "data" not in response or "hits" not in response["data"]:
            raise ValueError("No cases found in GDC response")
        
        cases = response["data"]["hits"]
        print(f"Found {len(cases)} cases")
        
        # Extract clinical information
        clinical_list = []
        for case in cases:
            case_id = case.get("id")
            demographic = case.get("demographic", {})
            diagnoses = case.get("diagnoses", [{}])[0] if case.get("diagnoses") else {}
            
            clinical_list.append({
                "case_id": case_id,
                "age_at_index": demographic.get("age_at_index"),
                "gender": demographic.get("gender"),
                "race": demographic.get("race"),
                "vital_status": demographic.get("vital_status"),
                "days_to_death": demographic.get("days_to_death"),
                "days_to_birth": demographic.get("days_to_birth"),
                "primary_diagnosis": diagnoses.get("primary_diagnosis"),
                "tumor_stage": diagnoses.get("tumor_stage"),
                "ajcc_pathologic_t": diagnoses.get("ajcc_pathologic_t"),
                "ajcc_pathologic_n": diagnoses.get("ajcc_pathologic_n"),
                "ajcc_pathologic_m": diagnoses.get("ajcc_pathologic_m"),
            })
        
        clinical_df = pd.DataFrame(clinical_list)
        clinical_path = self.data_dir / "clinical_data.tsv"
        clinical_df.to_csv(clinical_path, sep="\t", index=False)
        print(f"Saved clinical data to {clinical_path}")
        
        return clinical_df
    
    def get_pam50_subtypes(self, auto_download: bool = False) -> pd.DataFrame:
        # Get PAM50 subtype annotations for TCGA-BRCA cases, optionally auto-downloading if not found. Parameters: auto_download (bool) - If True, automatically downloads from TCGA publications URL if not found via API. Returns: pd.DataFrame - DataFrame with case_id and PAM50 subtype labels
        print("Querying GDC API for PAM50 subtype annotations...")
        
        # Method 1: Try annotations endpoint
        print("  Method 1: Checking annotations endpoint...")
        annotations_endpoint = "/annotations"
        params = {
            "filters": json.dumps({
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "cases.project.project_id",
                            "value": [self.PROJECT_ID]
                        }
                    },
                    {
                        "op": "in",
                        "content": {
                            "field": "annotation_type",
                            "value": ["PAM50"]
                        }
                    }
                ]
            }),
            "fields": "case_id,annotation_type,entity_id",
            "size": "10000",
            "format": "JSON"
        }
        
        try:
            response = self._query_gdc(annotations_endpoint, params)
            if "data" in response and "hits" in response["data"] and len(response["data"]["hits"]) > 0:
                annotations = response["data"]["hits"]
                print(f"  Found {len(annotations)} PAM50 annotations")
                annotation_list = []
                for ann in annotations:
                    annotation_list.append({
                        "case_id": ann.get("case_id"),
                        "annotation_type": ann.get("annotation_type"),
                        "entity_id": ann.get("entity_id")
                    })
                return pd.DataFrame(annotation_list)
        except Exception as e:
            print(f"  Annotations endpoint failed: {e}")
        
        # Method 2: Try to find supplementary files with PAM50 in name
        print("  Method 2: Searching for supplementary files with PAM50...")
        try:
            files_endpoint = "/files"
            params = {
                "filters": json.dumps({
                    "op": "and",
                    "content": [
                        {
                            "op": "in",
                            "content": {
                                "field": "cases.project.project_id",
                                "value": [self.PROJECT_ID]
                            }
                        },
                        {
                            "op": "in",
                            "content": {
                                "field": "files.data_category",
                                "value": ["Clinical"]
                            }
                        }
                    ]
                }),
                "fields": "file_id,file_name,file_size,cases.case_id",
                "size": "10000",
                "format": "JSON"
            }
            response = self._query_gdc(files_endpoint, params)
            if "data" in response and "hits" in response["data"]:
                files = response["data"]["hits"]
                pam50_files = [f for f in files if "pam50" in f.get("file_name", "").lower() or "subtype" in f.get("file_name", "").lower()]
                if pam50_files:
                    print(f"  Found {len(pam50_files)} potential PAM50 files")
                    file_list = []
                    for f in pam50_files:
                        file_list.append({
                            "file_id": f.get("id"),
                            "file_name": f.get("file_name"),
                            "case_id": f.get("cases", [{}])[0].get("case_id") if f.get("cases") else None
                        })
                    pam50_files_df = pd.DataFrame(file_list)
                    pam50_files_path = self.data_dir / "pam50_files_manifest.tsv"
                    pam50_files_df.to_csv(pam50_files_path, sep="\t", index=False)
                    print(f"  Saved PAM50 files manifest to {pam50_files_path}")
                    print("  Note: These files need to be downloaded and parsed manually")
        except Exception as e:
            print(f"  Supplementary files search failed: {e}")
        
        # Method 3: Check if PAM50 is in clinical data (sometimes it's there)
        print("  Method 3: Checking clinical data for subtype fields...")
        try:
            cases_endpoint = "/cases"
            params = {
                "filters": json.dumps({
                    "op": "in",
                    "content": {
                        "field": "cases.project.project_id",
                        "value": [self.PROJECT_ID]
                    }
                }),
                "fields": "case_id,diagnoses.molecular_subtype_method,diagnoses.morphology",
                "size": "100",
                "format": "JSON"
            }
            response = self._query_gdc(cases_endpoint, params)
            if "data" in response and "hits" in response["data"]:
                cases = response["data"]["hits"]
                # Check if any have subtype information
                has_subtype = False
                for case in cases[:5]:  # Check first 5 as sample
                    diagnoses = case.get("diagnoses", [{}])[0] if case.get("diagnoses") else {}
                    if diagnoses.get("molecular_subtype_method"):
                        has_subtype = True
                        break
                if has_subtype:
                    print("  Found potential subtype fields in clinical data")
                    print("  Note: Full clinical data may contain subtype information")
        except Exception as e:
            print(f"  Clinical subtype check failed: {e}")
        
        # If all methods fail, try auto-download if enabled
        if auto_download:
            print("\n  Method 4: Attempting automatic download from TCGA publications...")
            try:
                file_path = self.download_pam50_from_url()
                return self.load_pam50_from_file(str(file_path))
            except Exception as e:
                print(f"  Auto-download failed: {e}")
        
        # If all methods fail, provide guidance
        print("\n" + "=" * 60)
        print("PAM50 subtypes not found via GDC API")
        print("=" * 60)
        print("\nPAM50 subtypes are typically available from:")
        print("1. Supplementary files from TCGA publications:")
        print("   - BRCA.547.PAM50.SigClust.Subtypes.txt (from Nature 2012 paper)")
        print("   - Available at: https://tcga-data.nci.nih.gov/docs/publications/brca_2012/")
        print("   - Use loader.download_pam50_from_url() to download automatically")
        print("2. TCGAbiolinks R package:")
        print("   - Use TCGAquery_subtype(tumor='BRCA') function")
        print("   - Returns PAM50 labels for available samples")
        print("3. PanCancerAtlas data:")
        print("   - More complete PAM50 assignments")
        print("   - Available via TCGAbiolinks or direct download")
        print("\nFor Python implementation:")
        print("- Use loader.download_pam50_from_url() to download from TCGA publications")
        print("- Then use loader.load_pam50_from_file() to load the downloaded file")
        print("- Or use PAM50 gene signature to infer subtypes from expression data")
        print("=" * 60)
        
        return pd.DataFrame()
    
    def download_pam50_from_url(self, url: Optional[str] = None, file_name: str = "BRCA.547.PAM50.SigClust.Subtypes.txt") -> Path:
        # Download PAM50 subtypes file from TCGA publications URL. Parameters: url (str, optional) - URL to download from. If None, uses default TCGA publication URL, file_name (str) - Name to save the file as. Returns: Path - Path to downloaded file
        if url is None:
            # Try multiple URLs - the old NCI URL may redirect to HTML
            # Try CGHub/Synapse URL first (more reliable), then fallback
            # Try multiple possible URLs - the file location may vary
            alt_urls = [
                "https://gdc.cancer.gov/files/public/file/BRCA_PAM50_Subtypes.txt",  # GDC direct
                "https://gdc.cancer.gov/files/public/file/BRCA.547.PAM50.SigClust.Subtypes.txt",  # GDC with original name
                "https://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt",  # Original NCI
                "https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/BRCA.547.PAM50.SigClust.Subtypes.txt",  # Firehose
            ]
            url = alt_urls[0]
            alt_url = alt_urls[1] if len(alt_urls) > 1 else None
        else:
            alt_url = None
        
        file_path = self.data_dir / file_name
        
        # Check if file exists and is valid (not HTML)
        if file_path.exists():
            # Check if it's actually a data file, not HTML
            try:
                with open(file_path, 'rb') as f:
                    first_bytes = f.read(1024)
                    if first_bytes.startswith(b'<!DOCTYPE') or first_bytes.startswith(b'<html'):
                        print(f"Existing file appears to be HTML, will re-download...")
                    else:
                        print(f"PAM50 file already exists: {file_path}")
                        return file_path
            except:
                pass
        
        print(f"Downloading PAM50 subtypes from: {url}")
        print("Note: This file contains PAM50 subtype assignments for 547 TCGA-BRCA samples")
        
        try:
            response = self.session.get(url, stream=True, timeout=300, allow_redirects=True)
            response.raise_for_status()
            
            # Check content type to make sure it's not HTML
            content_type = response.headers.get('content-type', '').lower()
            if 'text/html' in content_type or 'html' in content_type:
                if alt_url:
                    print(f"  Received HTML instead of data file, trying alternative URL...")
                    response.close()
                    response = self.session.get(alt_url, stream=True, timeout=300, allow_redirects=True)
                    response.raise_for_status()
                else:
                    raise ValueError("Downloaded file appears to be HTML, not a data file")
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            with open(file_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            if downloaded % (1024 * 64) == 0:  # Print every 64KB
                                print(f"  Progress: {percent:.1f}%", end='\r')
            
            # Verify it's not HTML
            with open(file_path, 'rb') as f:
                first_bytes = f.read(1024)
                if first_bytes.startswith(b'<!DOCTYPE') or first_bytes.startswith(b'<html'):
                    raise ValueError("Downloaded file is HTML, not a data file. The URL may have changed.")
            
            print(f"\nDownloaded PAM50 subtypes to: {file_path}")
            return file_path
        except (requests.exceptions.RequestException, ValueError) as e:
            print(f"  Error downloading PAM50 file: {e}")
            print(f"\n  Alternative options:")
            print(f"  1. Manual download from TCGA publications:")
            print(f"     - Visit: https://gdc.cancer.gov/about-data/publications/brca_2012")
            print(f"     - Download: BRCA.547.PAM50.SigClust.Subtypes.txt")
            print(f"     - Save to: {self.data_dir / file_name}")
            print(f"  2. Use TCGAbiolinks R package to export PAM50 subtypes")
            print(f"  3. Check for updated URLs in TCGA documentation")
            raise
    
    def load_pam50_from_file(self, file_path: str) -> pd.DataFrame:
        # Load PAM50 subtypes from a local file (e.g., downloaded supplementary file). Parameters: file_path (str) - Path to PAM50 subtype file. Returns: pd.DataFrame - DataFrame with case_id/sample_id and PAM50 subtype
        print(f"Loading PAM50 subtypes from file: {file_path}")
        file_path_obj = Path(file_path)
        
        if not file_path_obj.exists():
            raise FileNotFoundError(f"PAM50 file not found: {file_path}")
        
        # Check if file is HTML (bad download)
        with open(file_path_obj, 'rb') as f:
            first_bytes = f.read(1024)
            if first_bytes.startswith(b'<!DOCTYPE') or first_bytes.startswith(b'<html'):
                print(f"  Error: File appears to be HTML, not a data file")
                print(f"  This usually means the download URL redirected to a webpage")
                print(f"  Please delete the file and try downloading again, or use an alternative source")
                return pd.DataFrame()
        
        # Try to parse common PAM50 file formats
        try:
            # Try TSV format first - skip bad lines
            df = pd.read_csv(file_path_obj, sep="\t", on_bad_lines='skip', comment='#', skip_blank_lines=True)
            
            # If empty or only one column, try with different separator
            if len(df.columns) == 1:
                df = pd.read_csv(file_path_obj, sep=None, engine='python', on_bad_lines='skip', comment='#', skip_blank_lines=True)
            
            # Remove completely empty rows
            df = df.dropna(how='all')
            
            print(f"  Loaded {len(df)} rows from file")
            print(f"  Columns: {list(df.columns)}")
            
            # Look for common column names
            subtype_cols = [col for col in df.columns if any(x in col.lower() for x in ['pam50', 'subtype', 'molecular', 'call'])]
            id_cols = [col for col in df.columns if any(x in col.lower() for x in ['case', 'sample', 'barcode', 'patient', 'tumor'])]
            
            # If no exact match, try first and second column
            if not subtype_cols or not id_cols:
                if len(df.columns) >= 2:
                    id_cols = [df.columns[0]]
                    subtype_cols = [df.columns[1]]
            
            if subtype_cols and id_cols:
                result = df[[id_cols[0], subtype_cols[0]]].copy()
                result.columns = ['case_id', 'pam50_subtype']
                # Clean up values
                result['case_id'] = result['case_id'].astype(str).str.strip()
                result['pam50_subtype'] = result['pam50_subtype'].astype(str).str.strip()
                result = result[result['pam50_subtype'].notna()]
                result = result[result['pam50_subtype'] != '']
                result = result[result['pam50_subtype'].str.lower() != 'nan']
                
                print(f"  Found {len(result)} samples with PAM50 subtypes")
                if len(result) > 0:
                    print(f"  Sample subtypes: {result['pam50_subtype'].value_counts().to_dict()}")
                return result
            else:
                print(f"  Warning: Could not identify subtype columns. Available columns: {list(df.columns)}")
                print(f"  First few rows:")
                print(df.head())
                return df
                
        except Exception as e:
            print(f"  Error parsing file: {e}")
            print(f"  File may be in an unexpected format")
            # Try to read first few lines to help debug
            try:
                with open(file_path_obj, 'r', encoding='utf-8', errors='ignore') as f:
                    lines = [f.readline().strip() for _ in range(10)]
                    print(f"  First few lines of file:")
                    for i, line in enumerate(lines[:5], 1):
                        print(f"    Line {i}: {line[:100]}")
            except:
                pass
            return pd.DataFrame()
    
    def download_file(self, file_id: str, file_name: str) -> Path:
        # Download a file from GDC by file ID. Parameters: file_id (str) - GDC file UUID, file_name (str) - Name to save the file as. Returns: Path - Path to downloaded file
        file_path = self.data_dir / file_name
        
        if file_path.exists():
            print(f"File already exists: {file_path}")
            return file_path
        
        print(f"Downloading {file_name}...")
        url = f"{self.GDC_API_BASE}/data/{file_id}"
        
        response = self.session.get(url, stream=True, timeout=300)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        if downloaded % (1024 * 1024) == 0:  # Print every MB
                            print(f"  Progress: {percent:.1f}%", end='\r')
        
        print(f"\nDownloaded: {file_path}")
        return file_path
    
    def load_expression_data(self, manifest_df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        # Load RNA-seq expression data. If files aren't downloaded, downloads them first. Parameters: manifest_df (pd.DataFrame, optional) - File manifest. If None, will query for it. Returns: pd.DataFrame - Expression matrix with genes as rows, samples as columns
        if manifest_df is None:
            manifest_df = self.get_file_manifest()
        
        # For now, return the manifest
        # Actual expression loading would require downloading and parsing files
        print("Note: Expression data loading requires downloading files.")
        print("Use download_file() to download specific files, then parse them.")
        
        return manifest_df
    
    def load_all(self, auto_download_pam50: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        # Load all TCGA-BRCA data: file manifest, clinical data, and PAM50 subtypes. Parameters: auto_download_pam50 (bool) - If True, automatically downloads PAM50 from TCGA publications if not found via API. Returns: tuple - (file_manifest, clinical_data, pam50_subtypes)
        print("=" * 60)
        print("Loading TCGA-BRCA data from GDC API")
        print("=" * 60)
        
        manifest = self.get_file_manifest()
        clinical = self.get_clinical_data()
        pam50 = self.get_pam50_subtypes(auto_download=auto_download_pam50)
        
        print("\n" + "=" * 60)
        print("Data loading complete!")
        print("=" * 60)
        
        return manifest, clinical, pam50
