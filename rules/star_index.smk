

GENOME_DIR = os.path.join(config["star_index_master_dir"], config["star_index_species"], config["star_index_version"], config["star_index_overhang_prefix"] + str(config["read_length"]), "")

if not os.path.exists(GENOME_DIR):
    os.system("mkdir -p {}".format(GENOME_DIR))

##

rule generate_genome:
    input:
        fasta = config["genome_fa"],
        gtf = config["gtf"]
    output:
        os.path.join(GENOME_DIR, "SA"),
        os.path.join(GENOME_DIR, "Genome")
    params:
        sjdbOverhang = config["read_length"] - 1
    threads:
        4
    conda:
        "../env/env_align.yaml"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {GENOME_DIR} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.sjdbOverhang} \
        --genomeSAsparseD 10
        """
