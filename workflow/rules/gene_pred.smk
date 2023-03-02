

rule select_mikado_train:
    input: "results/mikado/mikado.loci.metrics.tsv", "results/mikado/mikado.loci.gff3"
    output: "results/augustus_train/training.gff3"
    conda: "../envs/python3.yaml"
    params: f = 0.5, e = 2
    script: "../scripts/select_mik_train.py"


rule strict_repeat_filter_for_train:
    input: "results/augustus_train/training.gff3", config["genome"]+".out.gff"
    output: "results/augustus_train/training_filter.gff3"
    conda: "../envs/pybedtools.yaml"
    params: frac = 0.9
    script: "../scripts/filt_repeats_gtf2.py"


rule train_augustus:
    input: gff = "results/augustus_train/training_filter.gff3",
           g = config["genome"]+".masked"
    output: "results/augustus_training/autoAugTrain/training/test/augustus.2.CRF.out"
    conda: '../envs/augustus.yaml'
    params: spname = GENOME, odir = "results/augustus_training"
    shell: "mkdir -p {params.odir} && autoAugTrain.pl --trainingset={input.gff} --genome={input.g} --CRF "
           "--species={params.spname} --workingdir={params.odir} --optrounds=1 --verbose --useexisting"


rule portcullis_to_hint:
    input: "results/portcullis/portcullis.filtered.bam"
    output: "results/hints_for_augustus/intronhints.gff"
    conda: '../envs/augustus.yaml'
    shell: "bam2hints --intronsonly --minintronlen=15 --in={input} --out={output}"


rule mikado_to_hint:
    input: "results/mikado/mikado.loci.gtf"
    output: "results/hints_for_augustus/mikado_exons_hints.gff"   
    conda: "../envs/python3.yaml"
    params: t = "mikado", s = 'E', f = 'exonpart'
    script: "../scripts/gtf_to_hints.py"


if config['metaeuk_only']:
    rule metaeuk_full_annot:
        input: genome =  config["genome"]+'.masked', profile = config["prot_fasta"]
        output: out = "results/metaeuk.gff"
        params: output = lambda w, output: output[0].replace(".gff", ''), tmp = "results/metaeuk/tmp/"
        conda: "../envs/metaeuk.yaml"
        threads: 30
        shell: "mkdir -p {params.tmp} && metaeuk easy-predict {input.genome} {input.profile} {params.output} {params.tmp} --threads {threads} && "
               "rm -r {params.tmp}"
else:
    rule metaeuk:
        input: genome =  config["genome"]+'.masked', profile = config["prot_fasta"]
        output: out = "results/metaeuk/pred_results.gff"
        params: output = lambda w, output: output[0].replace(".gff", ''), tmp = "results/metaeuk/tmp/"
        conda: "../envs/metaeuk.yaml"
        threads: 30
        shell: "metaeuk easy-predict {input.genome} {input.profile} {params.output} {params.tmp} --threads {threads} && "
               "rm -r {params.tmp}"

rule metaeuk_to_gtf:
    input: "results/metaeuk/pred_results.gff"
    output: "results/metaeuk/pred_results.gtf"
    conda: "../envs/agat.yaml"
    shell: "agat_convert_sp_gxf2gxf.pl --gff {input} -o {output}"


rule metaeuk_to_hint:
    input: "results/metaeuk/pred_results.gtf"
    output: "results/hints_for_augustus/protsim_hints.gff"
    conda: "../envs/python3.yaml"
    params: t = "metaeuk", s = 'P', f = "CDSpart"
    script: "../scripts/gtf_to_hints.py"


rule cat_hints:
  input:
        "results/hints_for_augustus/intronhints.gff",
        "results/hints_for_augustus/mikado_exons_hints.gff", 
        "results/hints_for_augustus/protsim_hints.gff"
  output: "results/hints_for_augustus/hints.gff"
  shell: "cat {input} > {output}"


rule split_fasta:
    input: genome = config["genome"] + ".masked"
    output: temp(expand('.'.join(config["genome"].split('.')[:-1]) + ".{i}." + config["genome"].split('.')[-1] + ".masked",\
                        i=['0'+str(i) for i in range(0, 10)] + [str(i) for i in range(10, 20)]))
    conda: "../envs/pyfasta.yaml"
    shell: "pyfasta split {input} -n 20"


if not config['metaeuk_only']:
    rule augustus:
        input: hints = "results/hints_for_augustus/hints.gff",
               c = config.get("augustus_config"),
               train = "results/augustus_training/autoAugTrain/training/test/augustus.2.CRF.out",
               genome = '.'.join(config["genome"].split('.')[:-1]) + ".{i}." + config["genome"].split('.')[-1] + ".masked"
        output: f"results/augustus/{GENOME}.{{i}}.aug.out"
        params: sp = GENOME
        conda: '../envs/augustus.yaml'
        shell: "augustus --uniqueGeneId=true --gff3=on --species={params.sp} --hintsfile={input.hints} "
               "--extrinsicCfgFile={input.c} --allow_hinted_splicesites=atac "
               "--alternatives-from-evidence=false {input.genome} > {output}"


    rule cat_augustus:
        input: expand(f"results/augustus/{GENOME}.{{i}}.aug.out", i=['0'+str(i) for i in range(0, 10)] + [str(i) for i in range(10, 20)])
        output: f"results/augustus/{GENOME}.aug.gff3"
        shell: "cat {input} | grep -v '^#' > {output}"