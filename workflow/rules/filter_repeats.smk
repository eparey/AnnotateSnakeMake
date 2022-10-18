rule cds_to_fasta:
    input: gff = f"results/augustus/{GENOME}.aug.gff3",
           genome = config["genome"] + ".masked"
    output: "results/augustus/pep.fa"
    conda: "../envs/agat.yaml"
    shell: 'agat_sp_extract_sequences.pl --gff {input.gff} -f {input.genome} -o {output} -t CDS --merge -p'


#hmmpress -f {input.hmm} #hmm = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/Pfam-A.hmm"
rule pfam_scan:
    input: f = "results/augustus/pep.fa", db = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/Pfam-A.hmm.h3i"
    output: out = "results/pfam/unfiltered.genes.pfam-domains.txt"
    conda: "../envs/pfam.yaml"
    params: d = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/"
    threads: 8
    shell: "pfam_scan.pl -fasta {input.f} -dir {params.d} -outfile {output.out} -cpu {threads}"

rule busco:
    input: "results/augustus/pep.fa"
    output: "results/busco/busco_unfilterred/short_summary.specific.metazoa_odb10.busco_unfilterred.json"
    threads: 4
    conda: "../envs/busco.yaml"
    params: odir = "results/busco/"
    shell: "busco -l metazoa_odb10 --mode proteins --tar -o busco_unfilterred -f -i {input} --cpu {threads} --out_path {params.odir}"


rule filter_gene_models:
    input: f"results/augustus/{GENOME}.aug.gff3", f'results/repeats/{GENOME}_repeats.bed'
    output: "results/filter_models/gene_models.filt.{i}.gff3"
    conda: "../envs/pybedtools.yaml"
    params: frac = lambda w: w.i, cut = 0.1
    script: "../scripts/filt_repeats_gtf2.py"


rule plot_nbgenes:
    input: f = expand("results/filter_models/genelist.filt.{i}.txt",\
                      i=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    output: "results/filter_models/genes.svg"
    params: ','.join([f"filter_{i}" for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]])
    conda: "../envs/plots.yaml"
    script: "../scripts/nbgenes_plots.py"


rule cds_to_fasta_filtered:
    input: gff = "results/filter_models/gene_models.filt.{i}.gff3", genome = config["genome"] + ".masked"
    output: "results/filter_models/pep.filt{i}.fa"
    conda: "../envs/agat.yaml"
    shell: 'agat_sp_extract_sequences.pl --gff {input.gff} -f {input.genome} -o {output} -t CDS --merge -p'


rule pfam_scan_filt:
    input: f = "results/filter_models/pep.filt{i}.fa" , d = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/Pfam-A.hmm.h3i"
    output: "results/filter_models/pfam/filt{i}.genes.pfam-domains.txt"
    params: d = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/"
    conda: "../envs/pfam.yaml"
    shell: "pfam_scan.pl -fasta {input.f} -dir {params.d} -outfile {output} -cpu {threads}"


rule plot_pfam:
    input: nof = "results/pfam/unfiltered.genes.pfam-domains.txt",
           f = expand("results/filter_models/pfam/filt{i}.genes.pfam-domains.txt",
                      i=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    output: "results/filter_models/pfam.svg"
    conda: "../envs/plots.yaml"
    params: ','.join(["no_filter"] + [f"filter_{i}" for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]])
    script: "../scripts/pfam_plots.py"


rule busco_filt:
    input: "results/filter_models/pep.filt{i}.fa"
    output: "results/busco/busco_filter{i}/short_summary.specific.metazoa_odb10.busco_filter{i}.json"
    params: jname = lambda w: "busco_filter" + w.i, odir = "results/busco/"
    threads: 4
    conda: "../envs/busco.yaml"
    shell: "busco -l metazoa_odb10 --tar --mode proteins -o {params.jname} -f -i {input} --cpu {threads} --out_path {params.odir}"


rule plot_busco:
    input: nof = "results/busco/busco_unfilterred/short_summary.specific.metazoa_odb10.busco_unfilterred.json",
           f = expand("results/busco/busco_filter{i}/short_summary.specific.metazoa_odb10.busco_filter{i}.json",
                      i=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    output: "results/filter_models/busco.svg"
    conda: "../envs/plots.yaml"
    params:  ','.join(["no_filter"] + [f"filter_{i}" for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]])
    script: "../scripts/busco_plots.py"



# if config.get('repeats_filter', 'auto') == 'auto':
#     rule select_threshold:
#         input: 
#         output: "results/filter_models/gene_models.filtered.gff3"
#         script:

# else:
#     rule skip_select_threshold:
#         input: "results/filter_models/gene_models.filt.{i}.gff3"
#         output: "results/filter_models/gene_models.filtered.gff3"
#         shell: "cp {input} {output}"