=======
Single-Cell RNA-Seq Analysis Workflow (Kallisto/Bustools -> scVI)This repository contains a Snakemake workflow and associated scripts to process single-cell RNA-seq data starting from SRA files, using Kallisto/Bustools for quantification, and scvi-tools for integration and latent space generation. It also includes a template notebook for downstream analysis based on the scVI results.Workflow OverviewData Download: Downloads SRA data using prefetch.FASTQ Conversion: Converts SRA to FASTQ using fasterq-dump.Reference Indexing: Builds a Kallisto index from a provided transcriptome FASTA.Quantification: Runs kallisto bus to pseudoalign reads.Counting: Runs bustools (sort, correct, count) to generate gene x cell count matrices (MTX format).Preprocessing: Aggregates count matrices from all samples, performs QC filtering (min genes/cells, mitochondrial percentage), using Scanpy (scripts/preprocess_data_from_mtx.py).Annotation: Merges external metadata (cell types, conditions, etc.) from a CSV file with the aggregated data (scripts/annotate_adata.py).scVI Training: Trains an scVI model on the annotated, preprocessed counts, correcting for sample batch effects. Saves the trained model and the AnnData object with the learned latent space (scripts/run_scvi_training.py).Directory Structureyour_project_name/
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
Setup1. Clone Repository:# Replace <your-repo-url> with the actual URL
git clone <your-repo-url>
cd your_project_name
2. Create Conda Environment:Requires Conda (Miniconda or Anaconda) to be installed.conda env create -f envs/environment.yaml
conda activate snakemake_scvi
(Note: If using GPU with scVI, ensure your CUDA drivers are compatible and the correct PyTorch/CUDA versions are specified or installed in environment.yaml)3. Obtain Reference Files:Place the following files in the reference/ directory (or update paths in config.yaml):Transcriptome FASTA: e.g., Homo_sapiens.GRCh38.cdna.all.r112.fa.gz. Download from Ensembl or GENCODE. Make sure the release number matches the one in config.yaml.Transcript-to-Gene Map: e.g., transcripts_to_genes_ensembl112.txt. Generate using Ensembl BioMart (two columns: TranscriptIDGeneID, no header). Ensure it corresponds to the FASTA release.Barcode Whitelist: e.g., 3M-february-2018_v3.txt. Obtain the correct whitelist for your 10x chemistry (e.g., from a Cell Ranger download). Use the unzipped .txt version.4. Prepare Input Data Files:Place the following files in the data/ directory:SRR_Acc_List.txt: A text file with one SRA Run Accession ID per line for all your samples.ethnicity.csv (or rename): Your metadata file. Must contain a column with cell barcodes and columns for annotations (cell type, ethnicity, infection status, etc.). Ensure the column names match those specified in config/config.yaml under metadata_columns.5. Configure Workflow:CRITICALLY EDIT config/config.yaml:Update all paths for reference files (fasta, t2g_map, whitelist).Verify the kallisto_bus: technology: parameter (e.g., 10xv3).Verify the metadata_file path and the metadata_columns mapping (keys must match CSV headers, values are desired AnnData .obs names).Adjust threads, max_epochs, n_latent, QC parameters etc. as needed.Set scvi: use_gpu: to true or false.Running the Workflow with SnakemakeNote on Conda Errors: Users have reported issues running conda commands within Snakemake subprocesses on some systems ("conda: command not found"), even when the environment is activated. If you encounter this:Ensure conda init <your_shell> (e.g., conda init bash) has been run previously.Crucially, close and reopen your terminal after running conda init.Reactivate the environment (conda activate snakemake_scvi).If problems persist, consider the --conda-frontend mamba (if mamba is installed) or --conda-executable $(which conda) options for Snakemake, or running the steps manually (see below).1. Activate Environment:conda activate snakemake_scvi
2. Dry Run (Optional but Recommended):Check the plan without executing commands.snakemake -n --use-conda
3. Execute Workflow:Run the pipeline. Replace N with the number of CPU cores you want to use.snakemake --cores N --use-conda
Snakemake will execute the steps, downloading data, running alignment/counting, preprocessing, annotation, and training the scVI model. Logs for each step will be saved in the logs/ directory.Running Steps Manually (Alternative if Snakemake Fails)If Snakemake execution is problematic, you can run the core Python steps manually after completing the initial data download and Kallisto/Bustools quantification (Steps 1-5 in the Snakemake workflow, which could be run using shell scripts based on the PDF or the Snakemake rule commands). Ensure the output directories match the structure expected by the scripts (defined in config.yaml).Activate Environment:conda activate snakemake_scvi
Preprocessing (MTX -> QC'd H5AD):# This assumes bustools counts are in results/bustools_counts/SAMPLE_ID/
python scripts/preprocess_data_from_mtx.py -c config/config.yaml
Annotation (Merge Metadata):# This takes the output from step 2 and adds metadata from the CSV
python scripts/annotate_adata.py -c config/config.yaml
scVI Training:# This takes the output from step 3 and trains the model
python scripts/run_scvi_training.py -c config/config.yaml
(Check script arguments using -h if you need to override specific inputs/outputs defined in the config)OutputsThe primary final outputs (located in results/ by default, within subdirectories defined in config.yaml) are:results/scvi_model/adata_scvi_trained.h5ad: AnnData object containing the scVI latent space in .obsm['X_scVI'].results/scvi_model/model/: Directory containing the saved scVI model state.Intermediate files (FASTQ, BUS, MTX, preprocessed/annotated AnnData) will also be generated if run via Snakemake.Downstream AnalysisUse the output file results/scvi_model/adata_scvi_trained.h5ad for subsequent analysis.Open the notebooks/downstream_analysis_template.ipynb notebook. It demonstrates loading the final AnnData object and outlines the steps for clustering, UMAP,
>>>>>>> 3b5653d (Initial commit of Snakemake workflow structure)
