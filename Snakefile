# Snakefile
import pandas as pd
import os

# --- Load Configuration ---
configfile: "config/config.yaml"

# --- Get Sample IDs ---
try:
    with open(config["sample_file"], 'r') as f:
        SAMPLES = [line.strip() for line in f if line.strip()]
    if not SAMPLES:
        raise ValueError("Sample file is empty.")
except FileNotFoundError:
    raise FileNotFoundError(f"Sample file not found at: {config['sample_file']}")
except Exception as e:
    raise RuntimeError(f"Error reading sample file: {e}")

print(f"Found {len(SAMPLES)} samples: {SAMPLES[:5]}...")

# --- Target Rule (Define final outputs) ---
rule all:
    input:
        config["output"]["scvi_adata"],
        config["output"]["scvi_model_dir"] + "model.pt" # Check actual saved model filename

# --- Workflow Rules ---

# Rule to download SRA files (runs once for all samples)
# Note: Prefetch downloads to a central cache, not directly to output dir
rule download_sra:
    output:
        touch("results/sra_download_complete.flag") # Flag file to mark completion
    log:
        "logs/prefetch.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Prefetching SRA files...' && "
        "prefetch --option-file {config[sample_file]} &> {log} && "
        "echo 'Prefetch complete.'"

# Rule to convert SRA to FASTQ (per sample)
rule dump_fastq:
    input:
        sra_flag = rules.download_sra.output
    output:
        r1 = config["workdir"] + "/fastq/{sample}_1.fastq.gz",
        r2 = config["workdir"] + "/fastq/{sample}_2.fastq.gz"
    params:
        outdir = config["workdir"] + "/fastq",
        threads = config["fasterq_dump"]["threads"],
        split = "--split-files" if config["fasterq_dump"]["split_files"] else ""
    log:
        "logs/fasterq_dump/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Converting {wildcards.sample} to FASTQ...' && "
        "mkdir -p {params.outdir} && "
        "fasterq-dump {wildcards.sample} "
        "{params.split} "
        "--outdir {params.outdir} "
        "-e {params.threads} "
        "--gzip " # Compress output
        "&> {log} && "
        "echo 'FASTQ conversion complete for {wildcards.sample}'"

# Rule to build Kallisto index (runs once)
rule kallisto_index:
    input:
        fasta = config["reference"]["fasta"]
    output:
        index = config["reference"]["index"]
    log:
        "logs/kallisto_index.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Building Kallisto index...' && "
        "kallisto index -i {output.index} {input.fasta} &> {log} && "
        "echo 'Kallisto index built.'"

# Rule to run Kallisto BUS (per sample)
rule kallisto_bus:
    input:
        index = rules.kallisto_index.output.index,
        r1 = config["workdir"] + "/fastq/{sample}_1.fastq.gz",
        r2 = config["workdir"] + "/fastq/{sample}_2.fastq.gz"
    output:
        bus = config["workdir"] + "/kallisto_bus/{sample}/output.bus",
        ec = config["workdir"] + "/kallisto_bus/{sample}/matrix.ec",
        tsv = config["workdir"] + "/kallisto_bus/{sample}/run_info.json",
        tx_names = config["workdir"] + "/kallisto_bus/{sample}/transcripts.txt" # Needed by bustools count
    params:
        outdir = config["workdir"] + "/kallisto_bus/{sample}",
        tech = config["kallisto_bus"]["technology"],
        threads = config["kallisto_bus"]["threads"]
    log:
        "logs/kallisto_bus/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Running Kallisto BUS for {wildcards.sample}...' && "
        "mkdir -p {params.outdir} && "
        "kallisto bus -i {input.index} -o {params.outdir} -x {params.tech} -t {params.threads} {input.r1} {input.r2} &> {log} && "
        "echo 'Kallisto BUS complete for {wildcards.sample}'"

# Rule to sort BUS file (per sample)
rule bustools_sort:
    input:
        bus = rules.kallisto_bus.output.bus
    output:
        sorted_bus = config["workdir"] + "/bustools_sort/{sample}/output.sorted.bus"
    params:
        outdir = config["workdir"] + "/bustools_sort/{sample}",
        threads = config["bustools"]["threads"],
        mem = "4G" # Adjust memory if needed
    log:
        "logs/bustools_sort/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Sorting BUS file for {wildcards.sample}...' && "
        "mkdir -p {params.outdir} && "
        "bustools sort -t {params.threads} -m {params.mem} -o {output.sorted_bus} {input.bus} &> {log} && "
        "echo 'BUS sort complete for {wildcards.sample}'"

# Rule to correct barcodes (per sample)
rule bustools_correct:
    input:
        sorted_bus = rules.bustools_sort.output.sorted_bus,
        whitelist = config["bustools"]["whitelist"]
    output:
        corrected_bus = config["workdir"] + "/bustools_correct/{sample}/output.corrected.bus"
    params:
        outdir = config["workdir"] + "/bustools_correct/{sample}"
    log:
        "logs/bustools_correct/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Correcting barcodes for {wildcards.sample}...' && "
        "mkdir -p {params.outdir} && "
        "bustools correct -w {input.whitelist} -o {output.corrected_bus} {input.sorted_bus} &> {log} && "
        "echo 'Barcode correction complete for {wildcards.sample}'"

# Rule to generate count matrix (per sample)
rule bustools_count:
    input:
        corrected_bus = rules.bustools_correct.output.corrected_bus,
        ec = config["workdir"] + "/kallisto_bus/{sample}/matrix.ec",
        tx_names = config["workdir"] + "/kallisto_bus/{sample}/transcripts.txt",
        t2g = config["reference"]["t2g_map"]
    output:
        mtx = config["workdir"] + "/bustools_counts/{sample}/genes.mtx",
        barcodes = config["workdir"] + "/bustools_counts/{sample}/genes.barcodes.txt",
        genes = config["workdir"] + "/bustools_counts/{sample}/genes.genes.txt"
    params:
        outdir = config["workdir"] + "/bustools_counts/{sample}",
        prefix = config["workdir"] + "/bustools_counts/{sample}/genes" # Prefix for output files
    log:
        "logs/bustools_count/{sample}.log"
    conda:
        "envs/environment.yaml"
    shell:
        "echo 'Generating counts for {wildcards.sample}...' && "
        "mkdir -p {params.outdir} && "
        "bustools count --genecounts -o {params.prefix} -g {input.t2g} -e {input.ec} -t {input.tx_names} {input.corrected_bus} &> {log} && "
        "echo 'Count generation complete for {wildcards.sample}'"

# Rule to aggregate, preprocess, and QC all samples
rule aggregate_preprocess:
    input:
        # Expand requires list of samples; get MTX file for each sample
        mtx_files = expand(config["workdir"] + "/bustools_counts/{sample}/genes.mtx", sample=SAMPLES),
        sample_file = config["sample_file"] # Pass sample file to script
    output:
        h5ad = config["output"]["preprocessed_adata"]
    params:
        bustools_base_dir = config["workdir"] + "/bustools_counts",
        min_genes = config["qc"]["min_genes_per_cell"],
        min_cells = config["qc"]["min_cells_per_gene"],
        mito_cutoff = config["qc"]["mito_cutoff_percent"],
        mito_prefix = config["qc"]["mito_prefix"]
    log:
        "logs/aggregate_preprocess.log"
    conda:
        "envs/environment.yaml"
    script:
        "scripts/preprocess_data_from_mtx.py"

# Rule to annotate the preprocessed data with external metadata
rule annotate_adata:
    input:
        adata_in = config["output"]["preprocessed_adata"],
        metadata_csv = config["metadata_file"]
    output:
        adata_out = config["output"]["annotated_adata"]
    params:
        # Pass column mapping from config to the script
        # This requires the script to accept argparse arguments for these
        barcode_col = config["metadata_columns"]["Barcode"],
        celltype_col = config["metadata_columns"]["predicted.celltype.12"],
        ethnicity_col = config["metadata_columns"]["Ethnicity"],
        infection_col = config["metadata_columns"]["Infection"]
        # Add others if needed
    log:
        "logs/annotate_adata.log"
    conda:
        "envs/environment.yaml"
    script:
        "scripts/annotate_adata.py"


# Rule to train scVI model
rule scvi_train:
    input:
        # Use the *annotated* adata as input
        h5ad = config["output"]["annotated_adata"]
    output:
        h5ad_out = config["output"]["scvi_adata"],
        model_dir = directory(config["output"]["scvi_model_dir"]), # Mark as directory output
        model_flag = config["output"]["scvi_model_dir"] + "model.pt" # Check actual model file name
    params:
        analysis_name = os.path.basename(config["output"]["scvi_model_dir"].rstrip('/')), # e.g., "model"
        max_epochs = config["scvi"]["max_epochs"],
        batch_key = config["scvi"]["batch_key"],
        n_latent = config["scvi"]["n_latent"],
        use_gpu = config["scvi"]["use_gpu"] # Pass GPU preference
    log:
        "logs/scvi_train.log"
    conda:
        "envs/environment.yaml"
    script:
        "scripts/run_scvi_training.py"

