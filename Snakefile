import os
import snakemake.io
import glob
import sys


def prepare_rnaSpades(param, reads):
    # print(f"original: {reads}", file=sys.stderr)
    # print(f"str: {str(reads)}", file=sys.stderr)
    reads = list(reads)
    # print(f"list: {reads}", file=sys.stderr)
    all_params = list()
    for read in reads:
        all_params.append(f"{param} {read}")
    # print(' '.join(all_params), file=sys.stderr)
    return ' '.join(all_params)

ROOT_DIR = "/home/mabuelanin/dib-dev/denovo-rnaseq/workflow/"

SAMPLES_DIR = ROOT_DIR + "samples"
TRIMMED_SAMPLES = ROOT_DIR + "trimmed"
ASSEMBLY_DIR = ROOT_DIR + "assembly"
SALMON_QUANT = ROOT_DIR + "salmon_quant"
DESEQ2_OUT_DIR = ROOT_DIR + "DESEQ2"
DIFF_EXP_RESULTS = DESEQ2_OUT_DIR + "/complete"


SAMPLES, = glob_wildcards(SAMPLES_DIR + "/{sample}_1.fastq.gz")

rule all:
    input:
        # fastp trimming, this can be turned off in favor of previous or next rules
        expand("{OUTDIR}" + "/trimmed_{sample}_{file_type}.fastq.gz", OUTDIR = TRIMMED_SAMPLES, sample=SAMPLES, file_type = ["R1_PE", "R2_PE", "merged"]),
        
        # denovo assembly
        ASSEMBLY_DIR + "/transcripts.fasta",

        # salmon index
        SALMON_QUANT + "/transcripts.index/info.json",

        # salmon quantification
        expand(SALMON_QUANT + "/{sample}_quant/quant.sf", sample = SAMPLES),

        # aggregated quantification
        SALMON_QUANT + "/agg_quant/" + "aggr_quant.isoform.TMM.EXPR.matrix",

        # diff exp
        DIFF_EXP_RESULTS + "/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.count_matrix",

        # filter diff expressed genes P0.005_C1
        DIFF_EXP_RESULTS + "/_done_P0.005_C1"

rule extract_diff_expressed:
    input:
        count_matrix = DIFF_EXP_RESULTS + "/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.count_matrix",
        agg_quant = SALMON_QUANT + "/agg_quant/aggr_quant.isoform.TMM.EXPR.matrix",
        samples_list = ROOT_DIR + "samples.tsv",

    
    params:
        p_val = 0.005,
        log_fold_change = 1,
        deseq_dir = DIFF_EXP_RESULTS,

    output:
        DIFF_EXP_RESULTS + "/_done_P0.005_C1"

    shell: """
        OUT_DIR={params.deseq_dir}/FILTERED_P{params.p_val}-C{params.log_fold_change} && \
        mkdir -p $OUT_DIR && \
        cd {params.deseq_dir} && \
        TRINITY_HOME=$(ls -d -1 -tra $CONDA_PREFIX/../../pkgs/trinity-*/opt/trinity-* | tail -n 1) && \
        $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
        --matrix {input.agg_quant} \
        --samples {input.samples_list} \
        -P {params.p_val} -C {params.log_fold_change} && \
        mv *P{params.p_val}_C{params.log_fold_change}* $OUT_DIR && \
        touch _done_P{params.p_val}_C{params.log_fold_change}
    """


rule deseq2:
    input: 
        agg_quant = SALMON_QUANT + "/agg_quant/aggr_quant.isoform.TMM.EXPR.matrix",
        samples_list = ROOT_DIR + "samples.tsv",
    
    output:
        count_matrix = DIFF_EXP_RESULTS + "/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.count_matrix"

    params:
        deseq2_out_dir = DIFF_EXP_RESULTS,

    shell: """
    sed -i 's/_quant//g' {input.agg_quant} && \
    mkdir -p {params.deseq2_out_dir} && \
    TRINITY_HOME=$(ls -d -1 -tra $CONDA_PREFIX/../../pkgs/trinity-*/opt/trinity-* | tail -n 1) && \
    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
        --matrix {input.agg_quant} \
        --samples_file {input.samples_list} \
        --method DESeq2 \
        --output {params.deseq2_out_dir}
    """

rule aggregate_quants:
    input: 
        _salmon_quant = expand(SALMON_QUANT + "/{sample}_quant/quant.sf", sample = SAMPLES),
    
    output:
        agg_file = SALMON_QUANT + "/agg_quant/aggr_quant.isoform.TMM.EXPR.matrix"

    params:
        agg_quant_dir = SALMON_QUANT + "/agg_quant",
        quant_list = SALMON_QUANT + "/agg_quant" + "/quant_files.list",
        out_prefix = SALMON_QUANT + "/agg_quant/aggr_quant",

    shell: """
        mkdir -p {params.agg_quant_dir} && \
        TRINITY_HOME=$(ls -d -1 -tra $CONDA_PREFIX/../../pkgs/trinity-*/opt/trinity-* | tail -n 1) && \
        find $(pwd)/workflow/*_quant -name "quant.sf" | tee {params.quant_list} && \
        $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon \
        --out_prefix {params.out_prefix} \
        --name_sample_by_basedir \
        --quant_files {params.quant_list} \
        --gene_trans_map none
    """

rule quantification:
    threads: 16
    input:
        _salmon_index = SALMON_QUANT + "/transcripts.index/info.json",
        r1_pe = SAMPLES_DIR + "/{sample}_1.fastq.gz",
        r2_pe = SAMPLES_DIR + "/{sample}_2.fastq.gz",
    
    output:
        SALMON_QUANT + "/{sample}_quant/quant.sf",
    
    params:
        salmon_index = SALMON_QUANT + "/transcripts.index",
        salmon_quant_dir = SALMON_QUANT,
    
    shell: """
        salmon quant -i {params.salmon_index} -p 2 -l IU -1 <(gunzip -c {input.r1_pe}) -2 <(gunzip -c {input.r2_pe}) -o {params.salmon_quant_dir}/{wildcards.sample}_quant
    """

rule salmon_index:
    input: 
        transcripts = ASSEMBLY_DIR + "/transcripts.fasta",

    output: SALMON_QUANT + "/transcripts.index/info.json"
    params:
        salmon_dir = SALMON_QUANT + "/transcripts.index"

    shell: """
        mkdir -p {params.salmon_dir} && \
        salmon index --index {params.salmon_dir} --transcripts {input.transcripts}
    """



rule denovo_assembly:
    threads: 32

    input:
        R1_PE = expand(TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq.gz", sample = SAMPLES),
        R2_PE = expand(TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq.gz", sample = SAMPLES),
        MERGED = expand(TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz", sample = SAMPLES),

    output:
        transcripts = ASSEMBLY_DIR + "/transcripts.fasta",    
    params:
        spades_output_dir = ASSEMBLY_DIR,
        cores = 16,
        rnaspades_tmp_dir = ASSEMBLY_DIR + "/tmp",
        prepared_R1_PE = prepare_rnaSpades("-1", expand(TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq.gz", sample = SAMPLES)),
        prepared_R2_PE = prepare_rnaSpades("-2", expand(TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq.gz", sample = SAMPLES)),
        prepared_MERGED = prepare_rnaSpades("--merged", expand(TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz", sample = SAMPLES)),
    
    resources:
        mem_mb = 32000,
        nodes = 1,
        time = 2000,
        partition = "med2"
    
    shell: """
        /usr/bin/time -v \
        rnaspades.py {params.prepared_R1_PE} {params.prepared_R2_PE} {params.prepared_MERGED} -t {params.cores} \
        --tmp-dir {params.rnaspades_tmp_dir} \
        -o {params.spades_output_dir}
    """

rule fastp:
    threads: 32
    resources: 
        mem_mb=30000, 
        time_min=60, 
    input:
        r1_pe = SAMPLES_DIR + "/{sample}_1.fastq.gz",
        r2_pe = SAMPLES_DIR + "/{sample}_2.fastq.gz",

    output:
        OP_R1_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_SE.fastq.gz",
        OP_R2_SE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_SE.fastq.gz",
        OP_R1_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq.gz",
        OP_R2_PE=TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq.gz",
        OP_FAILED=TRIMMED_SAMPLES + "/trimmed_{sample}_failed.fastq.gz",
        OP_MERGED=TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz",

    shell: """
        fastp --in1 {input.r1_pe} --in2 {input.r2_pe} \
        --cut_front --cut_right --cut_window_size 4 \
        --cut_mean_quality 5 --trim_poly_x --length_required 25 \
        --low_complexity_filter --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --merge \
        --merged_out {output.OP_MERGED} --out1 {output.OP_R1_PE} \
        --out2 {output.OP_R2_PE} --unpaired1 {output.OP_R1_SE} \
        --unpaired2 {output.OP_R2_SE} --failed_out {output.OP_FAILED}
    """