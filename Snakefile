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


SIGS_OUTDIR = ROOT_DIR + "sigs"
cDBG_OUTDIR = ROOT_DIR + "cDBGk75_samples"
cDBG_ASSEMBLY_DIR = ROOT_DIR + "assembled_cDBGk75_samples"
TRIMMED_ASSEMBLY_DIR = ROOT_DIR + "assembled_trimmed_samples"
ITS_DB_DIR = ROOT_DIR + "its_db"
SPLITTED_FASTA_DIR = ROOT_DIR + "splitted_fasta"
BLASTN_RESULTS_DIR = ROOT_DIR + "blastn_results"
SEQCLEAN_TRIMMED_ASSEMBLY_DIR = ROOT_DIR + "seqclean_assembled_trimmed_samples"
BUSCO_DATASET_DIR = ROOT_DIR + "BUSCO_DATASET"
BUSCO_REPORTS = ROOT_DIR + "BUSCO_REPORTS"


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


        # # Generate cDBG for all trimmed samples {R1_PE, R2_PE, Merged}
        # expand("{OUTDIR}" + "/cDBG_k75_all_samples.{EXT}", OUTDIR = cDBG_OUTDIR, EXT = ["histo", "unitigs.fa"]),
        # # Merge all signatures
        # SIGS_OUTDIR + "/all_samples_k31.sig",
        # # Assemble the cDBG of all trimmed samples
        # cDBG_ASSEMBLY_DIR + "/transcripts.fasta",
        # TRIMMED_ASSEMBLY_DIR + "/transcripts.fasta",
        # # Download and create blastdb of the ITS database
        # # expand(ITS_DB_DIR + "/its_8.3.{EXT}", EXT = ['nhr', 'nin', 'nsq']),
        # ITS_DB_DIR + "/its_8.3.fa.ndb",
        # # cDBG blastn query on ITS
        # BLASTN_RESULTS_DIR + "/its_cDBG_all_samples.blastn",
        # # BLASTN_RESULTS_DIR + "/its_assembled_trimmed_samples.blastn",
        # # seqclean of assembled trimmed samples output
        # TRIMMED_ASSEMBLY_DIR + "/cleaned_transcripts.fasta",

        # BUSCO_REPORTS + "/CLEANED_TRANSCRIPT_TRIMMED/short_summary.specific.eukaryota_odb10.CLEANED_TRANSCRIPT_TRIMMED.txt",

        # # BLASTn assembled trimmed samples
        # # expand(BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/its_cleaned_assembled_trimmed_samples{SPLIT_PART}.blastn", SPLIT_PART = ["%03d" % (num,) for num in range(1,33,1)])
        # BLASTN_RESULTS_DIR + "/cleaned_assembled_trimmed" + "/merged_its_cleaned_assembled_trimmed_samples.blastn",
        # TRIMMED_ASSEMBLY_DIR + "/busco_config.ini"

rule aggregate_quants:
    input: 
        _salmon_quant = expand(SALMON_QUANT + "/{sample}_quant/quant.sf", sample = SAMPLES),
        samples_list = ROOT_DIR + "samples.tsv",
    
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