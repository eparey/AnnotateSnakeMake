from pathlib import Path


GENOME_PATH = config["genome"]
GENOME = Path(GENOME_PATH).stem

def get_genome_size(input_file):
    with open(input_file, 'r') as infile:
        for line in infile:
            if "total length:" in line:
                size = int(line.split()[2])
                return size