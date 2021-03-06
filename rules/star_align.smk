
# make sure the output folder for STAR exists before running anything
def return_parsed_extra_params(extra_params):
    # starting blank
    cmd = ""
    # for key in extra parameters
    for key in extra_params:
        # append the key value pair if it's a parameter that needed something
        if extra_params[key]:
            cmd += " --{0} {1}".format(key,extra_params[key])
        else:  # otherwise if it's parameter that's just a flag, append just the flag
            cmd += " --{0}".format(key)

    return cmd


star_outdir = os.path.join(OUTPUT_DIR, config["bam_outdir_name"])
multiqc_outdir = os.path.join(OUTPUT_DIR, config["multiqc_outdir_name"])

if not os.path.exists(star_outdir):
    os.system("mkdir -p {0}".format(star_outdir))

if not os.path.exists(multiqc_outdir):
    os.system("mkdir -p {0}".format(multiqc_outdir))


rule run_star_se:
    input:
        generated_index = os.path.join(GENOME_DIR, "SA"),
        generated_genome = os.path.join(GENOME_DIR, "Genome"),
        one = lambda wildcards: get_fastq(wildcards.sample, OPTIONS, OUTPUT_DIR)
    output:
        os.path.join(star_outdir, "{sample}.Aligned.sortedByCoord.out.bam"),
        os.path.join(star_outdir, "{sample}.SJ.out.tab"),
        os.path.join(star_outdir, "{sample}.Log.final.out"),

    params:
        extra_star_parameters = return_parsed_extra_params(config['extra_star_parameters']),
        genomeDir = GENOME_DIR,
        outTmpDir = os.path.join(star_outdir, "{sample}_tmpdir"),
        outputPrefix = os.path.join(star_outdir, "{sample}.")

    threads:
        4

    conda:
        "../env/env_align.yaml"

    # group:
    #     "align"

    shell:
        """
        STAR --genomeDir {params.genomeDir} \
        --readFilesIn {input.one} \
        --outFileNamePrefix {params.outputPrefix} \
        --readFilesCommand zcat --runThreadN {threads} \
        {params.extra_star_parameters} \
        --outTmpDir {params.outTmpDir}
        """


rule index_bams:
    input:
        os.path.join(star_outdir, "{sample}.Aligned.sortedByCoord.out.bam")

    output:
        os.path.join(star_outdir, "{sample}.Aligned.sortedByCoord.out.bam.bai")

    conda:
        "../env/env_align.yaml"

    # group:
    #     "align"

    shell:
        """
        samtools index {input}
        """

rule multiqc_star:
    input:
        expand(star_outdir + "{sample}.Aligned.sortedByCoord.out.bam", sample = SAMPLES)

    output:
        os.path.join(multiqc_outdir, "multiqc_report.html")

    params:
        dir = star_outdir,
        outdir = multiqc_outdir

    conda:
        "../env/env_align.yaml"

    shell:
        """
        multiqc -p -o {params.outdir} {params.dir}
        """
