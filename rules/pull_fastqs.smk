### This goes from bam to fastq (if needed)

def get_bam_to_pull(sample, options):
    '''
    Return path to bam file (from sample table) that want to extract fastqs
    '''

    if options[sample]["file_type"] == "bam" and options[sample]["realign"] == 1:
        return options[sample]["path"]

    else:
        return ""

def get_samtools_mate_flag(sample, options):
    '''
    Which flag should be passed to samtools view to extract all 1st/2nd mates of paired-end fragments?
    '''

    if options[sample]["file_type"] == "bam" and options[sample]["realign"] == 1:

        if options[sample]["which_mate"] == "r1":
            return 40

        elif options[sample]["which_mate"] == "r2":
            return 80

        else:
            raise ValueError("Invalid string for which_mate with sample {} - must be one of 'r1' or 'r2'".format(sample))
    else:
        return None

## This is for debugging purposes... Uncomment if having troubles
# for s in SAMPLES:
#     print("output of get_bam_to_pull for {}\n".format(s))
#     print("{}\n".format(get_bam_to_pull(s, OPTIONS)))
#


rule bam_to_fastq:
    input:
        lambda wildcards: get_bam_to_pull(wildcards.sample, OPTIONS)

    output:
        os.path.join(OUTPUT_DIR, config["fastq_outdir_name"], "{sample}.pulled.fastq.gz")

    params:
        which_mate = lambda wildcards: get_samtools_mate_flag(wildcards.sample, OPTIONS),
        temp_bam = os.path.join(OUTPUT_DIR, config["fastq_outdir_name"],"{sample}.bam")
        # Seems like samtools fastq doesn't accept input from STDIN
    conda:
        "../env/env_align.yaml"
    shell:
        """
        samtools view -h -f {params.which_mate} -U {input} > {params.temp_bam}
        samtools fastq -o {output} {params.temp_bam}
        rm {params.temp_bam}
        """
