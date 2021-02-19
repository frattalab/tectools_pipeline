### This goes from bam to fastq (if needed)

def get_bam_to_pull(sample, options):
    '''
    Return path to bam file (from sample table) that want to extract fastqs
    '''

    if options[sample]["file_type"] == "bam" and options[sample]["realign"] == 1:
        return options[sample]["path"]

    else:
        return ""


def join_fastqs(sample, options):
    '''
    Controls behaviour of concat_fastqs
    Returns 1 if fastqs should be concatenated (which_mate == "both" and strandedness == "unstranded")
    else returns 0
    '''

    if options[sample]["which_mate"] == "both" and options[sample]["strandedness"] == "unstranded":
        # Should concatenate all pulled fastqs
        return 1
    else:
        # Shouldn't concatenate pulled fastqs
        return 0


def fastqs_to_concat(sample, options, output_dir):
    '''
    Controls which files are joined together in  concat_fastqs
    If want sample to be realigned (file_type == bam & realign == 1)
        if which_mate set to "r1" - return path to pulled r1 file only
        if which_mate set to "r2" - return path to pulled r2 file only
        if which_mate set to "both" - return space separated string of all pulled fastqs

    Essentially - generates string passed to cat (which will generate 'pulled' fastq for STAR alignment)
    '''

    if options[sample]["file_type"] == "bam" and options[sample]["realign"] == 1:

        r1_target = os.path.join(output_dir, config["fastq_outdir_name"], sample + ".pulled.r1.fastq.gz")
        r2_target = os.path.join(output_dir, config["fastq_outdir_name"], sample + ".pulled.r2.fastq.gz")
        singletons_target = os.path.join(output_dir, config["fastq_outdir_name"], sample + ".singletons.fastq.gz")

        if options[sample]["which_mate"] == "r1":
            return r1_target

        elif options[sample]["which_mate"] == "r2":
            return r2_target

        elif options[sample]["which_mate"] == "both":

            # Sanity check to make sure combined files will also be considered unstranded (when eventually get to Tectool)
            if options[sample]["strandedness"] == "unstranded":
                return " ".join([r1_target, r2_target, singletons_target])
            else:
                raise ValueError("Sample {0} is incorrectly configured for alignment by combining all mates together - which_mate should be 'both' & strandedness 'unstranded'".format(sample))

        else:
            raise ValueError("Sample {0} has Invalid option for which_mate {1} - should be one of 'r1', 'r2' & 'both'".format(sample, options[sample]["which_mate"]))

    else:
        # No concatenation needed as not re-aligning after pulling from a BAM
        return ""


# def get_samtools_mate_flag(sample, options):
#     '''
#     Which flag should be passed to samtools view to extract all 1st/2nd mates of paired-end fragments?
#     '''
#
#     if options[sample]["file_type"] == "bam" and options[sample]["realign"] == 1:
#
#         if options[sample]["which_mate"] == "r1":
#             return 40
#
#         elif options[sample]["which_mate"] == "r2":
#             return 80
#
#         else:
#             raise ValueError("Invalid string for which_mate with sample {} - must be one of 'r1' or 'r2'".format(sample))
#     else:
#         return None

## This is for debugging purposes... Uncomment if having troubles
# for s in SAMPLES:
#     print("output of get_bam_to_pull for {}\n".format(s))
#     print("{}\n".format(get_bam_to_pull(s, OPTIONS)))
#


rule bam_to_fastq:
    input:
        lambda wildcards: get_bam_to_pull(wildcards.sample, OPTIONS)

    output:
        fq1 = temp(os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.pulled.r1.fastq.gz")),
        fq2 = temp(os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.pulled.r2.fastq.gz")),
        singletons = temp(os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.singletons.fastq.gz"))

    params:
        # which_mate = lambda wildcards: get_samtools_mate_flag(wildcards.sample, OPTIONS),
        temp_prefix = os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}"),

    conda:
        "../env/env_align.yaml"

    shell:
        """
        samtools sort -n -T {params.temp_prefix} {input} | \
        samtools fastq -1 {output.fq1} \
        -2 {output.fq2} \
        -s {output.singletons} -
        """


rule concat_fastqs:
    input:
        fq1 = os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.pulled.r1.fastq.gz"),
        fq2 = os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.pulled.r2.fastq.gz"),
        singletons = os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.singletons.fastq.gz")

    output:
        temp(os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.pulled.fastq.gz"))

    params:
        to_join = lambda wildcards: join_fastqs(wildcards.sample, OPTIONS),
        targets = lambda wildcards: fastqs_to_concat(wildcards.sample, OPTIONS, OUTPUT_DIR)

    shell:
        """
        if [ {params.to_join} == 1 ]; then
            # Want to cat fastqs together
            cat {params.targets} > {output}
        else
            # Just move target FQ to output file
            cp {params.target} {output}
        fi
        """
