import anndata as ad
import scanpy as sc # Often useful alongside anndata
import os
import pandas as pd

print(f"Anndata version: {ad.__version__}")
print(f"Scanpy version: {sc.__version__}")

# --- Files to Combine ---
# Using the paths you provided
file1 = "/l/users/nazira.dunbayeva/finalprojectcb703/results_kallisto_scvi/scvi_output/all_samples_preprocessed.h5ad"
file2 = "/l/users/nazira.dunbayeva/finalprojectcb703/results_kallisto_scvi/scvi_output/all_samples_preprocessed_2.h5ad"

# --- Define Output File ---
# <<< EDIT THIS PATH if you want to save the combined file somewhere else or with a different name
output_file = "/l/users/nazira.dunbayeva/finalprojectcb703/results_kallisto_scvi/scvi_output/all_samples_combined_preprocessed.h5ad"

print(f"Input file 1: {file1}")
print(f"Input file 2: {file2}")
print(f"Output file: {output_file}")

# --- Load AnnData objects ---
print("\nLoading AnnData objects...")
try:
    adata1 = ad.read_h5ad(file1)
    print(f"Loaded adata1 with shape: {adata1.shape}")
    adata2 = ad.read_h5ad(file2)
    print(f"Loaded adata2 with shape: {adata2.shape}")
except FileNotFoundError as e:
    print(f"\nError: Input file not found: {e}")
    print("Please ensure both H5AD files exist at the specified paths.")
    exit(1)
except Exception as e:
    print(f"\nError loading AnnData files: {e}")
    exit(1)

# --- Optional: Inspect before merging ---
# print("\nObs columns in adata1:", adata1.obs.columns)
# print("Obs columns in adata2:", adata2.obs.columns)
# print("Var names identical?", adata1.var.index.equals(adata2.var.index))
# It's good practice to ensure the .var (genes) DataFrames are identical if using join='outer'

# --- Concatenate AnnData objects ---
print("\nConcatenating AnnData objects...")
try:
    # axis=0: concatenate along observations (cells) - default
    # join='outer': keep genes present in *either* object (fill missing with 0/NaN) - default
    #     use join='inner' to keep only genes present in *both* objects.
    # label='origin_file': adds an 'origin_file' column to .obs indicating which file the cell came from (using '0' and '1')
    # index_unique='_': adds suffix (-1, -2) if any cell barcodes are duplicated between files
    adata_combined = ad.concat(
        [adata1, adata2],
        axis=0,
        join='outer',
        label='origin_file', # You can rename this label if desired
        index_unique='_'
    )
    print(f"Concatenation successful. Combined shape: {adata_combined.shape}")
    print(f"Added '.obs['origin_file']' to track source file ('0'='{os.path.basename(file1)}', '1'='{os.path.basename(file2)}')")

except Exception as e:
    print(f"\nError during concatenation: {e}")
    print("Concatenation failed. Check if .var indices (genes) are compatible or try join='inner'.")
    exit(1)

# --- Save Combined AnnData ---
print(f"\nSaving combined AnnData object to: {output_file}")
try:
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # Save the combined object
    adata_combined.write_h5ad(output_file, compression="gzip")
    print("Combined file saved successfully.")
except Exception as e:
    print(f"\nError saving combined AnnData file: {e}")
    exit(1)

print("\nScript finished.")