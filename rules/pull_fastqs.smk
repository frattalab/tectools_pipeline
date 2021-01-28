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



rule bam_to_fastq:
    input:
         get_bam_to_pull(wildcards.sample, OPTIONS)

    output:
        os.path.join(OUTPUT_DIR, fastq_outdir_name, "{sample}.pulled.fastq.gz")

    params:
        which_mate = get_samtools_mate_flag(wildcards.sample, OPTIONS)

    conda:
        "env/env_align.yaml"
    shell:
        """
        samtools view -f {params.which_mate} -U {input} | \
        samtools fastq -o {output} -
        """
