# scripts/preprocess_data_from_mtx.py

import scanpy as sc
import anndata as ad
import pandas as pd
import os
import argparse
import sys
import scipy.io
import traceback # Import added in previous debug step

print(f"Scanpy version: {sc.__version__}")
print(f"Anndata version: {ad.__version__}")
print(f"Pandas version: {pd.__version__}")

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="Load bustools MTX outputs, concatenate samples, calculate QC, filter, and save AnnData.")
parser.add_argument("-f", "--sample_file", required=True, help="Path to a file containing sample IDs, one per line.")
parser.add_argument("-b", "--bustools_base_dir", required=True, help="Base directory containing sample subdirectories with bustools count outputs (genes.mtx, etc.)")
parser.add_argument("-o", "--output_h5ad", required=True, help="Output path for the final preprocessed AnnData H5AD file.")
# --- New Optional Argument ---
parser.add_argument("-t", "--output_qc_txt", default="qc_counts_summary.txt", help="Output path for the QC cell counts summary text file (default: qc_counts_summary.txt)")
# ---------------------------
# Optional QC parameters
parser.add_argument("--min_genes", type=int, default=200, help="Minimum number of genes expressed per cell for filtering.")
parser.add_argument("--min_cells", type=int, default=3, help="Minimum number of cells expressing a gene for filtering.")
parser.add_argument("--mito_cutoff", type=float, default=15.0, help="Maximum mitochondrial percentage allowed per cell.")
parser.add_argument("--mito_prefix", type=str, default="MT-", help="Prefix for mitochondrial genes (e.g., 'MT-' for human, 'mt-' for mouse). Case-sensitive.")

args = parser.parse_args()

print("--- Parameters ---")
print(f"Sample File: {args.sample_file}")
print(f"Bustools Base Dir: {args.bustools_base_dir}")
print(f"Output H5AD: {args.output_h5ad}")
print(f"Output QC Txt: {args.output_qc_txt}") # Print new param
print(f"Min Genes/Cell: {args.min_genes}")
print(f"Min Cells/Gene: {args.min_cells}")
print(f"Mito Cutoff (%): {args.mito_cutoff}")
print(f"Mito Prefix: {args.mito_prefix}")
print("------------------")

# --- Read Sample IDs from File ---
# ... (Keep this section as it was in the previous working version) ...
print(f"Reading sample IDs from {args.sample_file}...")
try:
    with open(args.sample_file, 'r') as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    if not sample_ids:
        print(f"Error: No sample IDs found in {args.sample_file}. Exiting.")
        sys.exit(1)
    print(f"Found {len(sample_ids)} sample IDs: {', '.join(sample_ids[:10])}...")
except FileNotFoundError:
    print(f"Error: Sample file not found at {args.sample_file}. Exiting.")
    sys.exit(1)
except Exception as e:
     print(f"Error reading sample file {args.sample_file}: {e}. Exiting.")
     sys.exit(1)


# --- Data Loading Loop (Using explicit AnnData creation) ---
# ... (Keep the corrected loading loop from the previous working version) ...
adatas_list = []
print("\nLoading individual samples...")
for sample_id in sample_ids:
    print(f"  Loading {sample_id}...")
    counts_dir = os.path.join(args.bustools_base_dir, sample_id)
    mtx_file = os.path.join(counts_dir, "genes.mtx")
    barcodes_file = os.path.join(counts_dir, "genes.barcodes.txt")
    genes_file = os.path.join(counts_dir, "genes.genes.txt")

    if not all(os.path.exists(f) for f in [mtx_file, barcodes_file, genes_file]):
        print(f"  Warning: One or more input files (.mtx, .barcodes.txt, .genes.txt) not found for sample {sample_id} in {counts_dir}. Skipping.")
        continue

    try:
        matrix = scipy.io.mmread(mtx_file).tocsr()
        barcodes_df = pd.read_csv(barcodes_file, header=None, names=['barcode'])
        genes_df = pd.read_csv(genes_file, header=None, names=['gene_id'])
        obs_df = pd.DataFrame(index = sample_id + "_" + barcodes_df['barcode'].astype(str))
        obs_df['sample'] = sample_id
        var_df = pd.DataFrame(index = genes_df['gene_id'].astype(str))

        # Dimension check (keep from previous debug)
        if matrix.shape[0] == len(obs_df) and matrix.shape[1] == len(var_df):
             adata_sample = ad.AnnData(X=matrix, obs=obs_df, var=var_df) # No .T
             adata_sample.var_names_make_unique()
             print(f"  Loaded {sample_id} with shape {adata_sample.shape}")
             adatas_list.append(adata_sample)
        else:
             print(f"  Error: Dimension mismatch for {sample_id} before AnnData creation.")
             print(f"    Matrix shape C x G: {matrix.shape}")
             print(f"    Barcode list length C: {len(obs_df)}")
             print(f"    Gene list length G: {len(var_df)}")
             print(f"  Skipping {sample_id}.")
             continue

    except Exception as e:
        print(f"  Error processing sample {sample_id}: {e}. Skipping.")
        traceback.print_exc()

if not adatas_list:
    print("Error: No samples were successfully loaded into AnnData objects. Exiting.")
    sys.exit(1)


# --- Concatenate Samples Incrementally ---
# ... (Keep the incremental concatenation loop from the previous working version) ...
print(f"\nConcatenating {len(adatas_list)} samples incrementally...")
adata_combined = adatas_list[0].copy()
print(f"Started with sample {sample_ids[0]}, shape: {adata_combined.shape}")
for i in range(1, len(adatas_list)):
    current_sample_id = sample_ids[i]
    adata_to_add = adatas_list[i]
    print(f"Concatenating sample {current_sample_id} (index {i}) with shape {adata_to_add.shape}...")
    try:
        adata_combined = ad.concat(
            [adata_combined, adata_to_add], axis=0, label='sample',
            join='inner', index_unique=None, merge='unique'
        )
        print(f"  Successfully concatenated. New combined shape: {adata_combined.shape}")
    except Exception as e:
        print(f"  ERROR concatenating sample {current_sample_id}!")
        print(f"  Error details: {e}")
        traceback.print_exc()
        sys.exit(1)

print(f"\nFinal combined AnnData shape after concatenation: {adata_combined.shape}")

# --- Store cell count BEFORE QC ---
cells_before_qc = adata_combined.n_obs
print(f"\nTotal cells before QC filtering: {cells_before_qc}")
# --------------------------------

# --- QC Calculation ---
print("\nCalculating QC metrics...")
adata_combined.var['mt'] = adata_combined.var_names.str.startswith(args.mito_prefix)
sc.pp.calculate_qc_metrics(adata_combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print("QC metrics calculated.")

# --- QC Filtering ---
print("\nApplying QC filtering...")
initial_genes = adata_combined.n_vars # Store initial gene count for comparison

print(f"Filtering cells by min_genes = {args.min_genes}")
sc.pp.filter_cells(adata_combined, min_genes=args.min_genes)
print(f"Filtering genes by min_cells = {args.min_cells}")
sc.pp.filter_genes(adata_combined, min_cells=args.min_cells)

print(f"Filtering cells by mito_cutoff = {args.mito_cutoff}%")
mito_filter = adata_combined.obs.pct_counts_mt < args.mito_cutoff
print(f"  Cells passing mito filter: {mito_filter.sum()} / {adata_combined.n_obs}") # Note n_obs might have changed from gene filtering
adata_combined = adata_combined[mito_filter, :]

print(f"Finished filtering.")
print(f"  Shape after filtering: {adata_combined.shape}")

# --- Store cell count AFTER QC ---
cells_after_qc = adata_combined.n_obs
cells_removed = cells_before_qc - cells_after_qc
genes_after_qc = adata_combined.n_vars
genes_removed = initial_genes - genes_after_qc
print(f"Total cells remaining after QC filtering: {cells_after_qc}")
print(f"Total cells removed: {cells_removed}")
print(f"Total genes remaining after QC filtering: {genes_after_qc}")
print(f"Total genes removed: {genes_removed}")
# ---------------------------------

# --- Write QC counts to TXT file ---
if args.output_qc_txt:
    print(f"\nWriting QC summary to: {args.output_qc_txt}")
    try:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output_qc_txt), exist_ok=True)
        with open(args.output_qc_txt, 'w') as f:
            f.write(f"Cells before filtering: {cells_before_qc}\n")
            f.write(f"Cells after filtering: {cells_after_qc}\n")
            f.write(f"Cells removed by filtering: {cells_removed}\n")
            f.write(f"Genes before filtering: {initial_genes}\n")
            f.write(f"Genes after filtering: {genes_after_qc}\n")
            f.write(f"Genes removed by filtering: {genes_removed}\n")
        print("QC summary saved.")
    except Exception as e:
        print(f"  Warning: Could not write QC summary file: {e}")
# -----------------------------------


# --- Save Output AnnData ---
print(f"\nSaving final preprocessed AnnData object to: {args.output_h5ad}")
os.makedirs(os.path.dirname(args.output_h5ad), exist_ok=True)
adata_combined.write_h5ad(args.output_h5ad, compression="gzip")

print("Preprocessing script finished successfully.")