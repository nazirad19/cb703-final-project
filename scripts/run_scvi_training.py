# scripts/run_scvi_training.py

import scvi
import scanpy as sc
import anndata as ad
import argparse
import os
import torch # For checking GPU
import pandas as pd # For checking columns
import sys # For exit

print(f"SCVI version: {scvi.__version__}")
print(f"Scanpy version: {sc.__version__}")
print(f"Anndata version: {ad.__version__}")
print(f"Pandas version: {pd.__version__}")
print(f"PyTorch version: {torch.__version__}")

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="Train SCVI model on combined, preprocessed AnnData.")
parser.add_argument("-i", "--input_h5ad", required=True, help="Input preprocessed AnnData file (.h5ad)")
parser.add_argument("-o", "--output_h5ad", required=True, help="Output AnnData file with SCVI latent space (.h5ad)")
parser.add_argument("-m", "--output_model", required=True, help="Output directory to save the trained SCVI model")
parser.add_argument("-n", "--analysis_name", required=True, help="Name for this analysis run (used for logging)")
parser.add_argument("-e", "--max_epochs", type=int, default=400, help="Max training epochs for SCVI.")
parser.add_argument("-b", "--batch_key", type=str, default="sample", help="Key in adata.obs for batch information (e.g., 'sample')")
parser.add_argument("--n_latent", type=int, default=30, help="Number of latent dimensions for SCVI.")

args = parser.parse_args()

print("--- Parameters ---")
print(f"Input H5AD: {args.input_h5ad}")
print(f"Output H5AD: {args.output_h5ad}")
print(f"Output Model Dir: {args.output_model}")
print(f"Analysis Name: {args.analysis_name}")
print(f"Max Epochs: {args.max_epochs}")
print(f"Batch Key: {args.batch_key}")
print(f"N Latent: {args.n_latent}")
print("------------------")

# --- Load Data ---
print(f"Loading preprocessed data from {args.input_h5ad}...")
try:
    adata = sc.read_h5ad(args.input_h5ad)
    print(f"Loaded AnnData object with shape: {adata.shape}")
except FileNotFoundError:
    print(f"Error: Input file not found at {args.input_h5ad}")
    sys.exit(1)
except Exception as e:
    print(f"Error loading AnnData file: {e}")
    sys.exit(1)

# Ensure obs names are unique
adata.obs_names_make_unique()

# --- Check for necessary columns ---
if args.batch_key not in adata.obs.columns:
    print(f"Error: Specified batch key '{args.batch_key}' not found in adata.obs. Available columns: {adata.obs.columns.tolist()}")
    sys.exit(1)
# Check for columns needed for later interpretation/DE based on paper
for col in ['ethnicity', 'infection_status']:
     if col not in adata.obs.columns:
          print(f"Warning: Metadata column '{col}' (expected for downstream analysis based on paper) not found in adata.obs.")

# --- Check for GPU ---
use_gpu = torch.cuda.is_available()
print(f"Using GPU: {use_gpu}")

# --- Setup and Train SCVI ---
print("Setting up AnnData for SCVI...")
# Setup using sample ID as the batch key for technical correction.
# We do NOT include ethnicity or infection_status as covariates here,
# so we can analyze their effects later in the SCVI latent space.
# Specify raw counts layer if applicable (check your preprocessing script)
raw_counts_layer = None
if 'counts' in adata.layers:
    raw_counts_layer = 'counts'
    print(f"Using layer '{raw_counts_layer}' for SCVI training.")
else:
    print("Using adata.X for SCVI training (ensure these are raw/lightly normalized counts).")

scvi.model.SCVI.setup_anndata(
    adata,
    layer=raw_counts_layer,
    batch_key=args.batch_key
)

print("Initializing SCVI model...")
# Adjust n_latent, n_layers etc. if needed
model = scvi.model.SCVI(adata, n_latent=args.n_latent)
print(model)

print(f"Training SCVI model for up to {args.max_epochs} epochs...")
# Note: Removed use_gpu argument as per previous error resolution
model.train(max_epochs=args.max_epochs, early_stopping=True)
print("SCVI training complete.")

# --- Save Model ---
print(f"Saving trained model to: {args.output_model}")
# Ensure output directory exists
os.makedirs(args.output_model, exist_ok=True)
model.save(args.output_model, overwrite=True, save_anndata=False) # Don't save anndata within model dir

# --- Get Latent Space ---
print("Calculating latent representation...")
adata.obsm["X_scVI"] = model.get_latent_representation()
print("Added 'X_scVI' to adata.obsm")

# --- Save AnnData ---
print(f"Saving AnnData with latent space to: {args.output_h5ad}")
# Ensure output directory exists
os.makedirs(os.path.dirname(args.output_h5ad), exist_ok=True)
adata.write_h5ad(args.output_h5ad, compression="gzip")

print("--- SCVI Training Finished Successfully ---")