import os
import snakemake.io
import glob
import sys


def prepare_rnaSpades(param, reads):
    reads = list(reads)
    all_params = list()
    for read in reads:
        all_params.append(f"{param} {read}")
    return ' '.join(all_params)

ROOT_DIR = "/home/mhussien/hend/denovo-rnaseq/workflow/"

SAMPLES_DIR = ROOT_DIR + "samples"
TRIMMED_SAMPLES = ROOT_DIR + "trimmed"
ASSEMBLY_DIR = ROOT_DIR + "assembly"
SALMON_QUANT = ROOT_DIR + "salmon_quant"
DESEQ2_OUT_DIR = ROOT_DIR + "DESEQ2" 
TRANSDECODER_DIR = ROOT_DIR + "transdecoder"
TRINOTATE_DIR = ROOT_DIR + "trinotate"
TRINOTATE_DATA_DIR = ROOT_DIR + "trinotate/data"
TRINOTATE_OUT_DIR = ROOT_DIR + "trinotate/out_dir"


SAMPLES, = glob_wildcards(SAMPLES_DIR + "/{sample}_1.fastq.gz")

rule all:
    input:
        TRINOTATE_OUT_DIR + "/trinotate_annotation_report.xls",

        # # fastp trimming, this can be turned off in favor of previous or next rules
        # expand("{OUTDIR}" + "/trimmed_{sample}_{file_type}.fastq.gz", OUTDIR = TRIMMED_SAMPLES, sample=SAMPLES, file_type = ["R1_PE", "R2_PE", "merged"]),
        
        # # denovo assembly
        # ASSEMBLY_DIR + "/transcripts.fasta",

        # # salmon index
        # SALMON_QUANT + "/transcripts.index/info.json",

        # # salmon quantification
        # expand(SALMON_QUANT + "/{sample}_quant/quant.sf", sample = SAMPLES),

        # # aggregated quantification
        # SALMON_QUANT + "/agg_quant/" + "aggr_quant.isoform.TMM.EXPR.matrix",

        # # diff exp
        # DESEQ2_OUT_DIR + "/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.count_matrix",

        # # filter diff expressed genes P0.005_C1

        # # Extract expressed transcripts
        # DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/expressed_transcripts.fasta",

        # # Transdecoder
        # TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        # TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep",

        # # Trinotate Data
        # TRINOTATE_DATA_DIR + "/Pfam-A.hmm",
        # TRINOTATE_DIR + "/app_3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl",
        # TRINOTATE_DATA_DIR + "/uniprot_sprot.pep.pdb",
        # TRINOTATE_OUT_DIR + "/blastp.outfmt6",
        # TRINOTATE_OUT_DIR + "/blastx.outfmt6",
        # TRINOTATE_OUT_DIR + "/pfam.log",

        # TRINOTATE_OUT_DIR + "/gene_to_transcript.tsv",
        # TRINOTATE_OUT_DIR + "/trinotate_annotation_report.xls",



rule trinotate_report:
    threads: 1
    input:
        transdecoder_predicted = TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep",
        expressed_transcripts = DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/expressed_transcripts.fasta",
        sqlite = TRINOTATE_DATA_DIR + "/Trinotate.sqlite",
        blastp = TRINOTATE_OUT_DIR + "/blastp.outfmt6",
        blastx = TRINOTATE_OUT_DIR + "/blastx.outfmt6",
        trinotate_pfam = TRINOTATE_OUT_DIR + "/TrinotatePFAM.out",


    output:
        gene_to_transcript_map = TRINOTATE_OUT_DIR + "/gene_to_transcript.tsv",
        trinotate_report = TRINOTATE_OUT_DIR + "/trinotate_annotation_report.xls",

    params:
        trinotate_app = TRINOTATE_DIR + "/app_3.2.2/Trinotate",
        out_dir = TRINOTATE_OUT_DIR,

    resources:
        mem_mb = 25 * 1024,
        nodes = 1,
        time = 60 * 5,
        partition = "high2"

    shell: """
        set -e && \
        cat {input.expressed_transcripts} | grep ">" | cut -c2- |  awk -F_ '{{print $1 (NF>1? FS $2 : "")"\t"$0}}' > {output.gene_to_transcript_map} && \
        {params.trinotate_app} {input.sqlite} init --gene_trans_map {output.gene_to_transcript_map} --transcript_fasta {input.expressed_transcripts} --transdecoder_pep {input.transdecoder_predicted} && \
        {params.trinotate_app} {input.sqlite} LOAD_swissprot_blastp {input.blastp} && \
        {params.trinotate_app} {input.sqlite} LOAD_swissprot_blastx {input.blastx} && \
        {params.trinotate_app} {input.sqlite} LOAD_pfam {input.trinotate_pfam} && \
        {params.trinotate_app} {input.sqlite} report [opts] > {output.trinotate_report}
    """


rule identify_protein_domains:
    threads: 1
    input:
        blastp = TRINOTATE_OUT_DIR + "/blastp.outfmt6",
        blastx = TRINOTATE_OUT_DIR + "/blastx.outfmt6",
        transdecoder_predicted = TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep",
        pfam_A = TRINOTATE_DATA_DIR + "/Pfam-A.hmm"

    output:
        pfam_log = TRINOTATE_OUT_DIR + "/pfam.log",
        trinotate_pfam = TRINOTATE_OUT_DIR + "/TrinotatePFAM.out",
    
    params:
        out_dir = TRINOTATE_OUT_DIR
    
    resources:
        mem_mb = 25 * 1024,
        nodes = 1,
        time = 60 * 5,
        partition = "high2"

    shell: """
        cd {params.out_dir} && \
        hmmscan --cpu 12 \
        --domtblout TrinotatePFAM.out \
        {input.pfam_A} \
        {input.transdecoder_predicted} > {output.pfam_log}
    """

rule align_predicted_transdecoder:
    threads: 8
    input:
        _blast_db = TRINOTATE_DATA_DIR + "/uniprot_sprot.pep.pdb",
        transdecoder_predicted = TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep",
        
    
    params:
        blast_db = TRINOTATE_DATA_DIR + "/uniprot_sprot.pep"

    output:
        blastp = TRINOTATE_OUT_DIR + "/blastp.outfmt6",

    resources:
        mem_mb = 25 * 1024,
        nodes = 1,
        time = 60 * 5,
        partition = "high2"

    shell: """
       blastp -query {input.transdecoder_predicted} \
       -db {params.blast_db} \
       -num_threads 8 \
       -max_target_seqs 1 \
       -outfmt 6 \
       -evalue 1e-3 > {output.blastp}
    """

rule align_expressed_genes:
    threads: 8
    input:
        _blast_db = TRINOTATE_DATA_DIR + "/uniprot_sprot.pep.pdb",
        expressed_transcripts = DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/expressed_transcripts.fasta",
    
    params:
        blast_db = TRINOTATE_DATA_DIR + "/uniprot_sprot.pep"

    output:
        blastx = TRINOTATE_OUT_DIR + "/blastx.outfmt6",

    resources:
        mem_mb = 25 * 1024,
        nodes = 1,
        time = 60 * 5,
        partition = "high2"

    shell: """
        blastx -query {input.expressed_transcripts} \
        -db {params.blast_db} \
        -num_threads 8 \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-3 > {output.blastx}
    """
        


rule trinotate_prepare_blast:
    threads: 1
    input:
        uniprot = TRINOTATE_DATA_DIR + "/uniprot_sprot.pep",
    
    params:
        trinotate_data = TRINOTATE_DATA_DIR,

    output:
        TRINOTATE_DATA_DIR + "/uniprot_sprot.pep.pdb"
    
    resources:
        mem_mb = 20 * 1024,
        nodes = 1,
        time = 60 * 2,
        partition = "high2"
    
    shell: """
        makeblastdb -in {input.uniprot} -dbtype prot
    """


rule trinotate_prepare_DBs:
    threads: 1
    input:
        TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep",
        TRINOTATE_DIR + "/app_3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl",
    output:
        sqlite = TRINOTATE_DATA_DIR + "/Trinotate.sqlite",
        uniprot = TRINOTATE_DATA_DIR + "/uniprot_sprot.pep",
        pfam = TRINOTATE_DATA_DIR + "/Pfam-A.hmm",
    
    params:
        trinotate_dir = TRINOTATE_DIR,
        trinotate_data = TRINOTATE_DATA_DIR,
        trinotate_app = TRINOTATE_DIR + "/app_3.2.2"

    resources:
        mem_mb = 32 * 1024,
        nodes = 1,
        time = 5 * 60,
        partition = "high2"

    shell: """
        mkdir -p {params.trinotate_data} && \
        cd {params.trinotate_data} && \
        TRINOTATE_HOME=$CONDA_PREFIX/../../pkgs/trinotate*/bin/ && \
        {params.trinotate_app}/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate && \
        gunzip Pfam-A.hmm.gz && \
        hmmpress Pfam-A.hmm
    """

rule download_trinotate:
    threads: 1
    input:
        TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep",
    
    output:
        TRINOTATE_DIR + "/app_3.2.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl",

    params:
        trinotate_dir = TRINOTATE_DIR,
        trinotate_app = TRINOTATE_DIR + "/app_3.2.2"
    
    resources:
        mem_mb = 5 * 1024,
        nodes = 1,
        time = 5  * 60,
        partition = "high2"


    shell: """
    set -e
        mkdir -p {params.trinotate_dir} && \
        cd {params.trinotate_dir} && \
        wget https://github.com/Trinotate/Trinotate/archive/refs/tags/Trinotate-v3.2.2.zip && \
        unzip Trinotate-v3.2.2.zip && \
        echo mv Trinotate-Trinotate-v3.2.2 {params.trinotate_app} && \
        rm -rf {params.trinotate_app}
        mv Trinotate-Trinotate-v3.2.2 {params.trinotate_app} && \
        rm -rf tmp Trinotate-v3.2.2.zip
    """

rule detect_orf:
    threads: 1
    input:
        expressed_transcripts = DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/expressed_transcripts.fasta",

    params:
        transdecoder_dir = TRANSDECODER_DIR,

    output:
        TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        TRANSDECODER_DIR + "/expressed_transcripts.fasta.transdecoder.pep"

    resources:
        mem_mb = 10 * 1024,
        nodes = 1,
        time = 2 * 60,
        partition = "high2"

    shell: """
        mkdir -p {params.transdecoder_dir} && \
        cd {params.transdecoder_dir} && \
        TRANSDECODER_HOME=$(ls -d -1 -tra $CONDA_PREFIX/../../pkgs/transdecoder-*/opt/transdecoder | tail -n 1) && \
        $TRANSDECODER_HOME/TransDecoder.LongOrfs -t {input.expressed_transcripts} && \
        $TRANSDECODER_HOME/TransDecoder.Predict -t {input.expressed_transcripts}
    """


rule extract_expressed_genes:
    threads: 1
    input:
        assembled_transcripts = ASSEMBLY_DIR + "/transcripts.fasta",
        filtered_subset = DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.DE_results.P0.005_C1.DE.subset",
    
    output:
        out_fasta = DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/expressed_transcripts.fasta",
    
    resources:
        mem_mb = 1 * 1024,
        nodes = 1,
        time = 10,
        partition = "high2"

    shell: """
        seqkit grep -f \
        <(cat {input.filtered_subset} | sed '1d' | cut -f1) \
        {input.assembled_transcripts} -o {output.out_fasta}
    """



rule extract_diff_expressed:
    threads: 1
    input:
        count_matrix = DESEQ2_OUT_DIR + "/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.count_matrix",
        agg_quant = SALMON_QUANT + "/agg_quant/aggr_quant.isoform.TMM.EXPR.matrix",
        samples_list = ROOT_DIR + "samples.tsv",

    params:
        p_val = 0.005,
        log_fold_change = 1,
        deseq_dir = DESEQ2_OUT_DIR,

    output:
        filtered_subset = DESEQ2_OUT_DIR + "/FILTERED_P0.005-C1/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.DE_results.P0.005_C1.DE.subset",
    
    resources:
        mem_mb = 10 * 1024,
        nodes = 1,
        time = 60,
        partition = "high2"

    shell: """
        OUT_DIR={params.deseq_dir}/FILTERED_P{params.p_val}-C{params.log_fold_change} && \
        mkdir -p $OUT_DIR && \
        cd {params.deseq_dir} && \
        TRINITY_HOME=$(ls -d -1 -tra $CONDA_PREFIX/../../pkgs/trinity-*/opt/trinity-* | tail -n 1) && \
        $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
        --matrix {input.agg_quant} \
        --samples {input.samples_list} \
        -P {params.p_val} -C {params.log_fold_change} && \
        mv *P{params.p_val}_C{params.log_fold_change}* $OUT_DIR
    """


rule deseq2:
    threads: 1
    input: 
        agg_quant = SALMON_QUANT + "/agg_quant/aggr_quant.isoform.TMM.EXPR.matrix",
        samples_list = ROOT_DIR + "samples.tsv",
    
    output:
        count_matrix = DESEQ2_OUT_DIR + "/aggr_quant.isoform.TMM.EXPR.matrix.control_vs_treated.DESeq2.count_matrix"

    params:
        deseq2_out_dir = DESEQ2_OUT_DIR,
    
    resources:
        mem_mb = 10 * 1024,
        nodes = 1,
        time = 60,
        partition = "high2"

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
    threads: 1
    
    input: 
        _salmon_quant = expand(SALMON_QUANT + "/{sample}_quant/quant.sf", sample = SAMPLES),
    
    output:
        agg_file = SALMON_QUANT + "/agg_quant/aggr_quant.isoform.TMM.EXPR.matrix"

    params:
        agg_quant_dir = SALMON_QUANT + "/agg_quant",
        quant_list = SALMON_QUANT + "/agg_quant" + "/quant_files.list",
        out_prefix = SALMON_QUANT + "/agg_quant/aggr_quant",
    
    resources:
        mem_mb = 3 * 1024,
        nodes = 1,
        time = 60,
        partition = "high2"

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
    threads: 1

    input:
        _salmon_index = SALMON_QUANT + "/transcripts.index/info.json",
        r1_pe = SAMPLES_DIR + "/{sample}_1.fastq.gz",
        r2_pe = SAMPLES_DIR + "/{sample}_2.fastq.gz",
    
    output:
        SALMON_QUANT + "/{sample}_quant/quant.sf",
    
    params:
        salmon_index = SALMON_QUANT + "/transcripts.index",
        salmon_quant_dir = SALMON_QUANT,
    
    resources:
        mem_mb = 10 * 1024,
        nodes = 1,
        time = 60,
        partition = "high2"
    
    shell: """
        salmon quant -i {params.salmon_index} -p 2 -l IU -1 <(gunzip -c {input.r1_pe}) -2 <(gunzip -c {input.r2_pe}) -o {params.salmon_quant_dir}/{wildcards.sample}_quant
    """

rule salmon_index:
    threads: 1
    input: 
        transcripts = ASSEMBLY_DIR + "/transcripts.fasta",

    output: SALMON_QUANT + "/transcripts.index/info.json"
    params:
        salmon_dir = SALMON_QUANT + "/transcripts.index"
    
    resources:
        mem_mb = 10 * 1024,
        nodes = 1,
        time = 60,
        partition = "high2"

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
        cores = 32,
        rnaspades_tmp_dir = ASSEMBLY_DIR + "/tmp",
        prepared_R1_PE = prepare_rnaSpades("-1", expand(TRIMMED_SAMPLES + "/trimmed_{sample}_R1_PE.fastq.gz", sample = SAMPLES)),
        prepared_R2_PE = prepare_rnaSpades("-2", expand(TRIMMED_SAMPLES + "/trimmed_{sample}_R2_PE.fastq.gz", sample = SAMPLES)),
        prepared_MERGED = prepare_rnaSpades("--merged", expand(TRIMMED_SAMPLES + "/trimmed_{sample}_merged.fastq.gz", sample = SAMPLES)),
    
    resources:
        mem_mb = 100 * 1024,
        nodes = 1,
        time = 60 * 10,
        partition = "high2"
    
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