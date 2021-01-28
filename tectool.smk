import os
import pandas as pd

configfile: "config/config.yaml"


GTF = config["gtf"]

OUTPUT_DIR = os.path.join(config["master_output_dir"],"")

sample_tbl = pd.read_csv(config["sample_table_csv"], index_col="sample_name")

SAMPLES = sample_tbl.index.tolist()
OPTIONS = sample_tbl.to_dict(orient="index")

# print(sample_tbl)
# print(SAMPLES)
# print(OPTIONS)

##
include: "rules/pull_fastqs.smk"
include: "rules/star_index.smk"
include: "rules/star_index.smk"


rule all:
    input:
        expand(OUTPUT_DIR + "{sample}/enriched_annotation.gtf", sample = SAMPLES),
        expand(OUTPUT_DIR + "{sample}/classified_as_terminal_with_probabilities.tsv", sample = SAMPLES)

##


## These are helper functions to control which samples are sent for re-alignment, and how to find each sample's bam file


def get_bam(sample, options, output_dir):
    '''
    Returns path to input bam file for given sample
    If sample needed to be realigned, path will be output_dir/aligned/{sample}.bam
    If asked for sample to not be realigned, returns the path provided in the sample table/options

    params:
        sample <str>
        name of sample (in pipeline context should usually pass as wildcards.sample)

        options <dict>
        dict of {sample: {param1: x, param2: y}} generated from input sample table

        output_dir <str>
        path to master output directory (results for each sample stored within here)
    '''

    if not options[sample]["realign"] == 1 and options[sample]["file_type"] == "bam":

        return options[sample]["path"]

    else:

        return os.path.join(output_dir, "aligned", sample + ".bam")


def get_fastq(sample, options, output_dir):
    '''
    Return path to input fastq file for STAR aligning
    If sample has fastq file provided, return provided path
    If provided a bam file for realigning, this will return the path to the pulled fastq
    If provided a bam ready to go, this will return an empty string (don't need to run for this sample)

    params:
        sample <str>
        name of sample (in pipeline context should usually pass as wildcards.sample)

        options <dict>
        dict of {sample: {param1: x, param2: y}} generated from input sample table

        output_dir <str>
        path to master output directory (results for each sample stored within here)
    '''

    # if already provided fastq
    if options[sample][file_type] == "fastq":

        return options[sample][path]

    # Provided bam and need to pull out 1 of mates
    elif options[sample][file_type] == "bam" and options[sample][realign] == 1:

        return os.path.join(output_dir, "pulled_fastq", sample + ".fastq.gz")

    # Don't need to re-align
    else:
        return ""


rule run_tectool:
    input:
        bam = get_bam(wildcards.sample, OPTIONS, OUTPUT_DIR),
        bam_idx = get_bam(wildcards.sample, OPTIONS, OUTPUT_DIR) + ".bai",
        #blah blah blah

    output:
        os.path.join(OUTPUT_DIR, "{sample}/enriched_annotation.gtf"),
        os.path.join(OUTPUT_DIR, "{sample}/classified_as_terminal_with_probabilities.tsv")

    params:
        outdir = os.path.join(OUTPUT_DIR, "{sample}", ""),
        gtf = GTF,
        polya = config["polya_bed"],
        fa = config["genome_fa"],
        seq_dirn = OPTIONS[wildcards.sample]["strandedness"],
        min_jnc_reads = config["minimum_spliced_reads"],
        min_overlap = config["min_region_overlap"],
        max_fuzz = config["max_splice_fuzziness"],
        drop_selected = "--drop_manually_selected_features" if config["drop_manually_selected_features"] else "",
        drop_intronic = "--drop_intronic_polya_sites_of_overlapping_genes" if config["drop_intronic_overlap"] else ""

    conda:
        "env/env_tectool.yaml"

    shell:
        """
        tectool \
        --annotation {params.gtf} \
        --polyasites {params.polya} \
        --genome {params.fa} \
        --bam {input.bam} \
        --sequencing_direction {params.seq_dirn} \
        --minimum_spliced_reads_for_cryptic_exon_start_site {params.min_jnc_reads} \
        --min_region_overlap {params.min_overlap} \
        --max_splice_fuzziness {params.max_fuzz} \
        {params.drop_selected} \
        {params.drop_intronic} \
        --output_dir {params.output_dir} \
        --verbose
        """