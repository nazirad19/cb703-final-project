# Single-Cell RNA-seq Analysis of Ethnic Variations in Immune Response to falciparum Malaria

This repository contains a Snakemake workflow and associated scripts to process single-cell RNA-seq data starting from SRA files, using Kallisto/Bustools for quantification, and scvi-tools for integration and latent space generation. It also includes a template notebook for downstream analysis based on the scVI results.

## Authors

* Muhra Salem Mohamed Salem Almahri
* Esraa Saeed Abdalazeem Bayomi
* Nazira Dunbayeva

## Workflow Overview

1.  **Data Download:** Downloads SRA data using `prefetch`.
2.  **FASTQ Conversion:** Converts SRA to FASTQ using `fasterq-dump`.
3.  **Reference Indexing:** Builds a Kallisto index from a provided transcriptome FASTA.
4.  **Quantification:** Runs `kallisto bus` to pseudoalign reads.
5.  **Counting:** Runs `bustools` (sort, correct, count) to generate gene x cell count matrices (MTX format).
6.  **Preprocessing:** Aggregates count matrices from all samples, performs QC filtering (min genes/cells, mitochondrial percentage), using Scanpy (`scripts/preprocess_data_from_mtx.py`).
7.  **Annotation:** Merges external metadata (cell types, conditions, etc.) from a CSV file with the aggregated data (`scripts/annotate_adata.py`).
8.  **scVI Training:** Trains an scVI model on the annotated, preprocessed counts, correcting for sample batch effects. Saves the trained model and the AnnData object with the learned latent space (`scripts/run_scvi_training.py`).

## Directory Structure

```
your_project_name/
├── .gitignore
├── README.md
├── Snakefile
├── config/
│   └── config.yaml
├── data/
│   └── SRR_Acc_List.txt         # <-- User provides
│   └── ethnicity.csv            # <-- User provides
├── envs/
│   └── environment.yaml
├── notebooks/
│   └── downstream_analysis_template.ipynb
├── reference/                   # <-- User populates
│   ├── Homo_sapiens.GRCh38.cdna.all.rXXX.fa.gz
│   ├── transcripts_to_genes_ensemblXXX.txt
│   └── 10x_whitelist.txt
│   └── .gitignore
├── results/                     # <-- Created by Snakemake/scripts
│   └── .gitignore
├── logs/                        # <-- Created by Snakemake
│   └── .gitignore
└── scripts/
    ├── preprocess_data_from_mtx.py
    ├── annotate_adata.py
    └── run_scvi_training.py
```

## Setup

**1. Clone Repository:**
```bash
git clone https://github.com/nazirad19/cb703-final-project.git
cd cb703-final-project # Or your chosen directory name
```

**2. Create Conda Environment:**
Requires Conda (Miniconda or Anaconda) to be installed.
```bash
conda env create -f envs/environment.yaml
conda activate snakemake_scvi
```
*(Note: If using GPU with scVI, ensure your CUDA drivers are compatible and the correct PyTorch/CUDA versions are specified or installed in `environment.yaml`)*

**3. Obtain Reference Files:**
Place the following files in the `reference/` directory (or update paths in `config.yaml`):
* **Transcriptome FASTA:** e.g., `Homo_sapiens.GRCh38.cdna.all.r112.fa.gz`. Download from Ensembl or GENCODE. Make sure the release number matches the one in `config.yaml`.
* **Transcript-to-Gene Map:** e.g., `transcripts_to_genes_ensembl112.txt`. Generate using Ensembl BioMart (two columns: TranscriptID<TAB>GeneID, no header). Ensure it corresponds to the FASTA release.
* **Barcode Whitelist:** e.g., `3M-february-2018_v3.txt`. Obtain the correct whitelist for your 10x chemistry (e.g., from a Cell Ranger download). Use the unzipped `.txt` version.

**4. Prepare Input Data Files:**
Place the following files in the `data/` directory:
* `SRR_Acc_List.txt`: A text file with one SRA Run Accession ID per line for all your samples.
* `ethnicity.csv` (or rename): Your metadata file. Must contain a column with cell barcodes and columns for annotations (cell type, ethnicity, infection status, etc.). Ensure the column names match those specified in `config/config.yaml` under `metadata_columns`.

**5. Configure Workflow:**
* **CRITICALLY EDIT `config/config.yaml`**:
    * Update all paths for reference files (`fasta`, `t2g_map`, `whitelist`).
    * Verify the `kallisto_bus: technology:` parameter (e.g., `10xv3`).
    * Verify the `metadata_file` path and the `metadata_columns` mapping (keys must match CSV headers, values are desired AnnData `.obs` names).
    * Adjust `threads`, `max_epochs`, `n_latent`, QC parameters etc. as needed.
    * Set `scvi: use_gpu:` to `true` or `false`.

## Running the Workflow with Snakemake

**Note on Conda Errors:** Users have reported issues running `conda` commands within Snakemake subprocesses on some systems ("conda: command not found"), even when the environment is activated. If you encounter this:
1.  Ensure `conda init <your_shell>` (e.g., `conda init bash`) has been run previously.
2.  **Crucially, close and reopen your terminal** after running `conda init`.
3.  Reactivate the environment (`conda activate snakemake_scvi`).
4.  If problems persist, consider the `--conda-frontend mamba` (if mamba is installed) or `--conda-executable $(which conda)` options for Snakemake, or running the steps manually (see below).

**1. Activate Environment:**
```bash
conda activate snakemake_scvi
```

**2. Dry Run (Optional but Recommended):**
Check the planned job execution.
```bash
snakemake -n --use-conda
```

**3. Execute Workflow:**
Run the pipeline. Replace `N` with the number of CPU cores.
```bash
snakemake --cores N --use-conda
```

## Running Steps Manually (Alternative if Snakemake Fails)

If Snakemake execution is problematic, you can run the core Python steps manually after completing the initial data download and Kallisto/Bustools quantification (Steps 1-5 in the Snakemake workflow, which could be run using shell scripts based on the PDF or the Snakemake rule commands). Ensure the output directories match the structure expected by the scripts (defined in `config.yaml`).

1.  **Activate Environment:**
    ```bash
    conda activate snakemake_scvi
    ```
2.  **Preprocessing (MTX -> QC'd H5AD):**
    ```bash
    # Assumes bustools counts are in results/bustools_counts/SAMPLE_ID/
    # Adjust paths in config or script arguments if needed
    python scripts/preprocess_data_from_mtx.py -c config/config.yaml
    ```
3.  **Annotation (Merge Metadata):**
    ```bash
    # Takes the output from step 2 and adds metadata from the CSV
    python scripts/annotate_adata.py -c config/config.yaml
    ```
4.  **scVI Training:**
    ```bash
    # Takes the output from step 3 and trains the model
    python scripts/run_scvi_training.py -c config/config.yaml
    ```
    *(Check script arguments using `-h` if you need to override specific inputs/outputs defined in the config)*

## Outputs

The primary final outputs (located in `results/` by default, within subdirectories defined in `config.yaml`) are:

* `results/scvi_model/adata_scvi_trained.h5ad`: AnnData object containing the scVI latent space in `.obsm['X_scVI']`.
* `results/scvi_model/model/`: Directory containing the saved scVI model state.

Intermediate files (FASTQ, BUS, MTX, preprocessed/annotated AnnData) will also be generated if run via Snakemake.

## Downstream Analysis

Use the output file `results/scvi_model/adata_scvi_trained.h5ad` for subsequent analysis.

Open the `notebooks/downstream_analysis_template.ipynb` notebook. It demonstrates loading the final AnnData object and outlines the steps for clustering, UMAP, DE tests (using Scanpy), and generating the specific plots (proportions, volcano, heatmap, dot plot) based on the interactive session conducted previously. You will need to run the cells in the notebook and potentially adjust parameters or gene lists (especially for the dot plot) as needed.
