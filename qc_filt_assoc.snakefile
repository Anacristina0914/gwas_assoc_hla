import os
import pandas as pd
import shutil
onsuccess:
    shutil.rmtree(".snakemake")

# Inputs
# =====================================================
PHENO_FILE = config["pheno_file"]
PHENO_COLS = config["phenotypes"]
BFILE = config["bfile"]
OUTDIR = config["out_dir"]

# Output prefixes
# ======================================================
FILTERED_QC = os.path.join(OUTDIR, "01_filt", "qc_filt")
ASSOC   = os.path.join(OUTDIR, "02_assoc")
MODELS = os.path.join(OUTDIR, "03_models")
HLA_DIR = os.path.join(OUTDIR, "04_filt_assoc", "hla")
AA_DIR  = os.path.join(OUTDIR, "04_filt_assoc", "aa")
INDEL_DIR = os.path.join(OUTDIR, "04_filt_assoc", "indels")
PLOTS_DIR = os.path.join(OUTDIR, "05_assoc_plots", "aa")
SUMMARY_DIR = os.path.join(OUTDIR, "06_assoc_summary", "aa")
#LD    = os.path.join(OUTDIR, "03_LD_calc")

rule all:
     input:
        expand(os.path.join(ASSOC, "nocovar", "{pheno}.assoc.adjusted"), pheno=PHENO_COLS),
        expand(os.path.join(MODELS, "{pheno}.model"), pheno=PHENO_COLS),
        expand(os.path.join(ASSOC, "covar", "{pheno}.assoc.logistic"), pheno=PHENO_COLS),
        expand(os.path.join(HLA_DIR, "{pheno}.hla.assoc.logistic"), pheno=PHENO_COLS),
        expand(os.path.join(AA_DIR,  "{pheno}.aa.assoc.logistic"), pheno=PHENO_COLS),
        expand(os.path.join(INDEL_DIR, "{pheno}.indels.assoc.logistic"), pheno=PHENO_COLS),
        expand(os.path.join(AA_DIR, "{pheno}.aa_labs.assoc.logistic"), pheno=PHENO_COLS),
        expand(os.path.join(HLA_DIR, "{pheno}.hla_labs.assoc.logistic"), pheno=PHENO_COLS),
        expand(os.path.join(INDEL_DIR, "{pheno}.indels_labs.assoc.logistic"), pheno=PHENO_COLS),
        os.path.join(PLOTS_DIR, "assoc_typeI.jpeg"),
        os.path.join(PLOTS_DIR, "assoc_typeII.jpeg"),
        os.path.join(PLOTS_DIR, "top_hm_assoc_typeI.jpeg"),
        os.path.join(PLOTS_DIR, "top_hm_assoc_typeII.jpeg"),
        os.path.join(SUMMARY_DIR, "assoc_summary.tab")
        

# Step 1: QC filtering
# =========================================================
rule run_qc_filt:
    input:
        bed = BFILE + ".bed",
        bim = BFILE + ".bim",
        fam = BFILE + ".fam",
        filt_in = config["filt_in"]
    conda:
        "plink_env"
    output:
        bed = FILTERED_QC + ".bed",
        bim = FILTERED_QC + ".bim",
        fam = FILTERED_QC + ".fam"
    params:
        maf  = config["maf"],
        geno = config["geno"],
        hwe  = config["hwe"],
        filt = FILTERED_QC,
        bfile_prefix = BFILE
    log:
        os.path.join(OUTDIR, "logs", "plink_qc_filter.log")
    shell:
        """
        mkdir -p $(dirname {params.filt})
        mkdir -p $(dirname {log})
        plink --bfile {params.bfile_prefix} \
            --keep {input.filt_in} \
            --maf {params.maf} \
            --geno {params.geno} \
            --hwe {params.hwe} \
            --make-bed \
            --out {params.filt} > {log} 2>&1
        """
# Step 2: Associations
# ==========================================================
rule plink_assoc:
    input:
        bed = FILTERED_QC + ".bed",
        bim = FILTERED_QC + ".bim",
        fam = FILTERED_QC + ".fam",
        pheno = PHENO_FILE
    conda:
        "plink_env"
    output:
        assoc = os.path.join(ASSOC, "nocovar", "{pheno}.assoc"),
        adjust = os.path.join(ASSOC, "nocovar", "{pheno}.assoc.adjusted")
    params:
        ci = config["assoc_ci"],
        out = os.path.join(ASSOC, "nocovar", "{pheno}"),
        bfile_prefix = FILTERED_QC
    log:
        os.path.join(OUTDIR, "logs", "assoc_{pheno}.log")
    shell:
        """
        mkdir -p {ASSOC}/nocovar
        plink --bfile {params.bfile_prefix} \
            --pheno {input.pheno} \
            --pheno-name {wildcards.pheno} \
            --assoc --ci {params.ci} --adjust \
            --out {params.out} > {log} 2>&1
        """
rule plink_logist_assoc:
    input:
        bed = FILTERED_QC + ".bed",
        bim = FILTERED_QC + ".bim",
        fam = FILTERED_QC + ".fam",
        pheno = PHENO_FILE
    conda:
        "plink_env"
    output:
        assoc = os.path.join(ASSOC, "covar", "{pheno}.assoc.logistic"),
        adjust = os.path.join(ASSOC, "covar", "{pheno}.assoc.logistic.adjusted")
    params:
        ci = config["assoc_ci"],
        out = os.path.join(ASSOC, "covar", "{pheno}"),
        bfile_prefix = FILTERED_QC,
        adjust_vars = config["adjust_vars"]
    log:
        os.path.join(OUTDIR, "logs", "logist_assoc_{pheno}.log")
    shell:
        """
        mkdir -p {ASSOC}/covar
        plink --bfile {params.bfile_prefix} \
            --pheno {input.pheno} \
            --pheno-name {wildcards.pheno} \
            --logistic --ci {params.ci} --adjust \
            --covar {input.pheno} \
            --covar-name {params.adjust_vars} \
            --out {params.out} > {log} 2>&1
        """

# Step 3: Model test per phenotype
# ===========================================================
rule plink_model:
    input:
        bed = FILTERED_QC + ".bed",
        bim = FILTERED_QC + ".bim",
        fam = FILTERED_QC + ".fam",
        pheno = PHENO_FILE
    conda:
        "plink_env"
    output:
        model = os.path.join(MODELS, "{pheno}.model")
    params:
        out = os.path.join(MODELS, "{pheno}"),
        bfile_prefix = FILTERED_QC
    log:
        os.path.join(OUTDIR, "logs", "model_{pheno}.log")
    shell:
        """
        mkdir -p {MODELS}
        plink --bfile {params.bfile_prefix} \
              --pheno {input.pheno} \
              --pheno-name {wildcards.pheno} \
              --model \
              --out {params.out} > {log} 2>&1
        """

# Step 5.1: Extract only HLA alleles according to user-provided regex
# ==================================================================
rule extract_hla:
    input:
        assoc = os.path.join(ASSOC, "covar", "{pheno}.assoc.logistic")
    output:
        hla = os.path.join(HLA_DIR, "{pheno}.hla.assoc.logistic")
    conda:
        "plink_env"
    params:
        pattern = lambda wildcards: config["hla_pattern"],
        hla_dir = lambda wildcards: HLA_DIR
    log:
        os.path.join(OUTDIR, "logs", "hla_{pheno}.log")
    shell:
        """
        mkdir -p {params.hla_dir}
        mkdir -p $(dirname {log})
        awk 'NR==1 || $2 ~ /{params.pattern}/ && $5=="ADD"' {input.assoc} > {output.hla} 2> {log}
        echo "Extracted $(wc -l < {output.hla}) HLA variants for {params.pattern}" >> {log}
        """
# head -n 1 {input.assoc} > {output.hla}

# Step 5.2: Extract aminoacids according to user-provided regex
# =================================================================
rule extract_aa:
    input:
        assoc = os.path.join(ASSOC, "covar", "{pheno}.assoc.logistic")
    output:
        aa = os.path.join(AA_DIR, "{pheno}.aa.assoc.logistic")
    conda:
        "plink_env"
    params:
        pattern = lambda wildcards: config["aa_pattern"],
        aa_dir  = lambda wildcards: AA_DIR
    log:
        os.path.join(OUTDIR, "logs", "aa_{pheno}.log")
    shell:
        """
        mkdir -p {params.aa_dir}
        mkdir -p $(dirname {log})
        head -n 1 {input.assoc} > {output.aa}
        awk 'NR==1 || ($2 ~ /{params.pattern}/ && $2 !~ /_x$/ && $5=="ADD")' {input.assoc} > {output.aa} 2> {log}
        echo "Extracted $(wc -l < {output.aa}) AA variants for {params.pattern}" >> {log}
        """
## Step 5.3: Extract insertions and deletions from association files
# =================================================================
rule extract_indels:
    input:
        assoc = os.path.join(ASSOC, "covar", "{pheno}.assoc.logistic")
    output:
        indels = os.path.join(INDEL_DIR, "{pheno}.indels.assoc.logistic")
    conda:
        "plink_env"
    params:
        indel_dir = lambda wildcards: INDEL_DIR
    log:
        os.path.join(OUTDIR, "logs", "indels_{pheno}.log")
    shell:
        """
        mkdir -p {params.indel_dir}
        mkdir -p $(dirname {log})
        awk 'NR==1 || $2 ~ /INS/' {input.assoc} > {output.indels} 2> {log}
        awk '$2 ~ /x$/' {input.assoc} >> {output.indels} 2>> {log}
        echo "Extracted $(wc -l < {output.indels}) indel variants for {wildcards.pheno}" >> {log}
        """ 

# Step 6: Make allele labels for plotting.
# =====================================================================
rule make_labels:
    input:
        aa_filt = os.path.join(AA_DIR, "{pheno}.aa.assoc.logistic"),
        hla_filt = os.path.join(HLA_DIR, "{pheno}.hla.assoc.logistic"),
        indel_filt = os.path.join(INDEL_DIR, "{pheno}.indels.assoc.logistic")
    output:
        aa_labs = os.path.join(AA_DIR, "{pheno}.aa_labs.assoc.logistic"),
        hla_labs = os.path.join(HLA_DIR, "{pheno}.hla_labs.assoc.logistic"),
        indel_labs = os.path.join(INDEL_DIR, "{pheno}.indels_labs.assoc.logistic")
    conda:
        "plink_env"
    run:
        import pandas as pd
        # Make labels for aminoacids by parsing the SNP column
        aa_df = pd.read_csv(input.aa_filt, sep=r"\s+")
        aa_df[["type", "locus", "position", "genomic", "residue"]] = (
            aa_df["SNP"].str.split("_", expand=True, n=5)
            )
        aa_df["aa"] = aa_df.apply(
            lambda row: row["A1"] if pd.isna(row["residue"]) else row["residue"],
            axis=1
            )
        aa_df.to_csv(output.aa_labs, index=None, sep="\t")

        # Make labels for hla alleles by parsing the SNP column
        hla_df = pd.read_csv(input.hla_filt, sep=r"\s+")
        hla_df[["HLA", "gene", "allele"]] = (
            hla_df["SNP"].str.split("_", expand=True, n=2)
            )
        hla_df.to_csv(output.hla_labs, index=None, sep="\t") 

        # Make labels for indels by parsing the SNP column
        indel_df = pd.read_csv(input.indel_filt, sep=r"\s+")
        indel_df[["variant","gene", "position", "insertion_del"]] = (
            indel_df["SNP"].str.split("_", expand=True, n=3)
            )
        indel_df.to_csv(output.indel_labs, index=None, sep="\t") 


# Step 7: Make association plots for aminoacids
# =======================================================================
rule plot_aa_assoc:
    input:
        aa_labs = expand(os.path.join(AA_DIR, "{pheno}.aa_labs.assoc.logistic"), pheno=PHENO_COLS)
    output:
        type1_aa_plot = os.path.join(PLOTS_DIR, "assoc_typeI.jpeg"),
        type2_aa_plot = os.path.join(PLOTS_DIR, "assoc_typeII.jpeg"),
        hm_type1_plot = os.path.join(PLOTS_DIR, "top_hm_assoc_typeI.jpeg"),
        hm_type2_plot = os.path.join(PLOTS_DIR, "top_hm_assoc_typeII.jpeg"),
        summary_table = os.path.join(SUMMARY_DIR, "assoc_summary.tab")
    params:
        group_cols = config["group_cols"],
        top_n = config["top_n"]
    script:
        config["scripts"]["aa_assoc_plot"]