snakemake --configfile ../QuickAnnotMglacialis/config.yaml --directory ../QuickAnnotMglacialis/ --use-conda --cores 9 --until plot_busco

snakemake --config genome=Patella.vulgata_xgPatVulg1.1.fasta repeats_only=y --use-conda --cores 8 --until repeat_landscape --directory ~/projects/molluscs/genomes/Patella.vulgata/