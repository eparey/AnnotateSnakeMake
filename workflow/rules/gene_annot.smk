#FIXME: add a post-deploy in conda to run pasa 

rule format_aug_for_pasa:
    input: f"results/augustus/{GENOME}.aug.gff3"
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
    input: p = f"results/augustus/{GENOME}.aug.ok-for-pasa.gff3", g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config.txt"), d = "results/pasa/.created_db"
    output: touch("results/pasa/.loaded_augustus")
    conda: "../envs/pasa.yaml"
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/scripts/Load_Current_Gene_Annotations.dbi -c {input.c} -g {input.g} -P {input.p}"


checkpoint pasa_run_1:
    input: g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config_comp.txt"), t =  "results/mikado/mikado.subloci.fasta.clean", l = "results/pasa/.loaded_augustus"
    output: "GenomeMik.sqlite.gene_structures_post_PASA_updates.{i}.gff3"
    conda: "../envs/pasa.yaml"
    threads: 16
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c {input.c} -A -g {input.g} -t {input.t} --CPU {threads}"


rule pasa_update_db:
    input: p = checkpoints.pasa_run_1.output[0], g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config.txt")
    output: touch("results/pasa/.updated_db_after_run1")
    conda: "../envs/pasa.yaml"
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/scripts/Load_Current_Gene_Annotations.dbi -c {input.c} -g {input.g} -P {input.p}"

checkpoint pasa_run_2:
    input: up = "results/pasa/.updated_db_after_run1", g = config["genome"] + ".masked", c = config.get("pasa_conf", "../AnnotateSnakeMake/config/pasa_config_comp.txt"), t =  "results/mikado/mikado.subloci.fasta.clean",
    output: "GenomeMik.sqlite.gene_structures_post_PASA_updates.2760253.gff3"
    conda: "../envs/pasa.yaml"
    threads: 16
    shell: "${{CONDA_PREFIX}}/opt/pasa-2.5.2/Launch_PASA_pipeline.pl -c {input.c} -A -g {input.g} -t {input.t} --CPU {threads}"


rule clean_after_pasa:
    input: run1 =  checkpoints.pasa_run_1.output[0], run2 =  checkpoints.pasa_run_2.output[0]
    output: run1 = f"results/pasa/{GENOME}_pasa1_tmp.gff3", run2 = f"results/pasa/{GENOME}_pasa_gene_models_final.gff3"
    params: odir = lambda w, output: os.path.dirname(output[0])
    shell: "mv {input.run1} {output.run1} && mv {input.run2} {output.run2} && "
           "rm -r  __pasa_GenomeMik.sqlite* && mv GenomeMik.sqlite* {params.odir} "
           "mv alignment.validations.output {params.odir} && mv pasa_run.log.dir logs/"