if config['metaeuk_only']:
    rule format_metaeuk:
        input: "results/metaeuk.gff"
        output: "results/metaeuk.ok.gff"
        script: "../scripts/format_metaeuk.py"

    rule metaeuk_cds_to_fasta:
        input: gff = "results/metaeuk.ok.gff",
               genome = f"{GENOME_PATH}.masked"
        output: out = "results/augustus/pep.fa", tmp = "results/augustus/pep.tmp.fa" #dummy awk to uniquely identify seq
        conda: "../envs/pasa.yaml"
        shell: """
            gffread -y {output.tmp} -g {input.genome} {input.gff} && awk '{{if($0~"^>") {{print $0NR}} else{{print $0}}}}' {output.tmp} | sed s#/#-#g > {output.out}
        """
else:
    rule cds_to_fasta:
        input: gff = f"results/augustus/{GENOME}.aug.gff3",
               genome = f"{GENOME_PATH}.masked"
        output: "results/augustus/pep.fa"
        conda: "../envs/pasa.yaml"
        shell: """
            gffread -y {output} -g {input.genome} {input.gff}
        """


rule pfam_scan:
    input: f = "results/augustus/pep.fa", db = config['pfam_db'].rstrip('/') + '/Pfam-A.hmm.h3i'
    output: out = "results/pfam/unfiltered.genes.pfam-domains.txt"
    conda: "../envs/pfam.yaml"
    params: d =  config['pfam_db']
    threads: 8
    shell: "pfam_scan.pl -fasta {input.f} -dir {params.d} -outfile {output.out} -cpu {threads}"

rule busco:
    input: "results/augustus/pep.fa"
    output:
        out = "results/busco/busco_unfilterred/short_summary.specific.metazoa_odb10.busco_unfilterred.json"
    threads: 4
    conda: "../envs/busco.yaml"
    params: odir = "results/busco/"
    shell: "busco -l metazoa_odb10 --mode proteins --tar -o busco_unfilterred -f -i {input} --cpu {threads} --out_path {params.odir}"


if config['metaeuk_only']:


    rule filter_gene_models_meta:
        input: "results/metaeuk.ok.gff", f'results/repeats/{GENOME}_repeats.bed'
        output: "results/filter_models/gene_models.filt.{i}.gff3"
        conda: "../envs/pybedtools.yaml"
        params: frac = lambda w: w.i, cut = 0.1
        script: "../scripts/filt_repeats_gtf2.py"
else:
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
    input: gff = "results/filter_models/gene_models.filt.{i}.gff3", genome = f"{GENOME_PATH}.masked"
    output: out = "results/filter_models/pep.filt{i}.fa", tmp = "results/filter_models/pep.filt{i}.tmp.fa"
    conda: "../envs/pasa.yaml"
    shell: """
        gffread -y {output.tmp} -g {input.genome} {input.gff} && awk '{{if($0~"^>") {{print $0NR}} else{{print $0}}}}' {output.tmp} | sed s#/#-#g > {output.out}
    """

rule busco_filt:
    input: p = "results/filter_models/pep.filt{i}.fa",
           db = "results/busco/busco_unfilterred/short_summary.specific.metazoa_odb10.busco_unfilterred.json" #dummy to ensure db is downloaded
    output:
        out = "results/busco/busco_filter{i}/short_summary.specific.metazoa_odb10.busco_filter{i}.json"
    params: jname = lambda w: "busco_filter" + w.i, odir = "results/busco/"
    threads: 4
    conda: "../envs/busco.yaml"
    shell: "busco -l metazoa_odb10 --mode proteins --tar -o {params.jname} -f -i {input.p} --cpu {threads} --out_path {params.odir}"


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