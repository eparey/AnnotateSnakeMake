#FIXME: add a post-deploy in conda to run pasa without calling conda_prefix

rule format_aug_for_pasa:
    input:  "results/filter_models/gene_models.filt.0.6.gff"#f"results/augustus/gene_models_filter0.6.gff3" #f"results/augustus/{GENOME}.aug.gff3"
    output: f"results/augustus/{GENOME}.aug.ok-for-pasa.gff3"
    script: "../scripts/aug_to_pasa.py"

rule prepare_mikado_transcriptome_for_pasa:
    input: g = config["genome"] + ".masked", t = "results/mikado/mikado.subloci.gff3"
    output: fa = "results/mikado/mikado.subloci.fasta", cln = "results/mikado/mikado.subloci.fasta.clean"
    conda: "../envs/pasa.yaml"
    params: odir = lambda w, output: os.path.dirname(output[0])
    threads: 6
    shell: "gffread -w {output.fa} -g {input.g} {input.t} && "
           "${{CONDA_PREFIX}}/opt/pasa-2.5.2/bin/seqclean {output.fa} -c {threads} && "
           "mv mikado.subloci.* {params.odir} && rm -r cleaning_*/ &&"
           "mv err_seqcl_mikado.subloci.fasta.log logs/ && mv seqcl_mikado.subloci.fasta.log logs/"

rule pasa_prepare_db:
    input: g = config["genome"] + ".masked", t = "results/mikado/mikado.subloci.fasta.clean",
           u = "results/mikado/mikado.subloci.fasta", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config.txt")
    output: "results/pasa/GenomeMik.sqlite"
    conda: "../envs/pasa.yaml"
    threads: 16
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c {input.c} -C -r -R -g {input.g} -t {input.t} "
           "-T -u {input.u} --ALIGNERS blat,gmap --CPU {threads}"

rule pasa_load:
    input: p = f"results/augustus/{GENOME}.aug.ok-for-pasa.gff3", g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config.txt"), d = "results/pasa/GenomeMik.sqlite"
    output: touch("results/pasa/.loaded_augustus")
    conda: "../envs/pasa.yaml"
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/scripts/Load_Current_Gene_Annotations.dbi -c {input.c} -g {input.g} -P {input.p}"


rule pasa_run_1:
    input: g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config_comp.txt"), t =  "results/mikado/mikado.subloci.fasta.clean", l = "results/pasa/.loaded_augustus"
    output: f"results/pasa/pasa1_tmp.gff3"
    conda: "../envs/pasa.yaml"
    threads: 16
    shell: "rm GenomeMik.sqlite.gene_structures_post_PASA_updates.*.gff3 || true && "
           "${{CONDA_PREFIX}}/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c {input.c} -A -g {input.g} -t {input.t} --CPU {threads} && "
           "mv GenomeMik.sqlite.gene_structures_post_PASA_updates.*.gff3 {output}"

rule pasa_update_db:
    input: p = f"results/pasa/pasa1_tmp.gff3", g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config.txt")
    output: touch("results/pasa/.updated_db_after_run1")
    conda: "../envs/pasa.yaml"
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/scripts/Load_Current_Gene_Annotations.dbi -c {input.c} -g {input.g} -P {input.p}"

rule pasa_run_2:
    input: up = "results/pasa/.updated_db_after_run1", g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config_comp.txt"), t =  "results/mikado/mikado.subloci.fasta.clean"
    output: f"results/pasa/pasa_gene_models.gff3"
    conda: "../envs/pasa.yaml"
    threads: 16
    shell: "rm GenomeMik.sqlite.gene_structures_post_PASA_updates.*.gff3 || true && "
           "${{CONDA_PREFIX}}/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c {input.c} -A -g {input.g} -t {input.t} --CPU {threads} && "
           "mv GenomeMik.sqlite.gene_structures_post_PASA_updates.*.gff3 {output}"

rule clean_after_pasa:
    input: run2 =  "results/pasa/pasa_gene_models.gff3"
    output: touch("results/pasa/.end")
    # params: odir = lambda w, output: os.path.dirname(output[0]) + '/tmp/'
    shell: "rm -r  __pasa_GenomeMik.sqlite* || true && rm GenomeMik.sqlite* || true && "
           "rm alignment.validations.output || true && mv pasa_run.log.dir logs/ || true && "
           "rm blat.spliced_alignments.gff3 && rm 11.ooc && rm tmp-* && rm -r pblat_outdir/ &&"
           "rm gmap.spliced_alignments.gff* "