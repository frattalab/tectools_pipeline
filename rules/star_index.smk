

# star_index_master_dir: /SAN/vyplab/vyplab_reference_genomes/STAR/
# star_index_species: human/
# star_index_version: ensembl_GRCh38/
# star_index_overhang_prefix: star_indices_overhang
# read_length: 75

GENOME_DIR = os.path.join(config["star_index_master_dir"], config["star_index_species"], config["star_index_version"], config["star_index_overhang_prefix"] + str(config["read_length"]), "")

if not os.path.exists(GENOME_DIR):
    os.system("mkdir -p {}".format(GENOME_DIR))

##
