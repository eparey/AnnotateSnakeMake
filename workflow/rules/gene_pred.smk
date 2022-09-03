
#rule exclude_repeats

# rule prepare_taco_output:
#     input: "results/transdecoder/transcripts.fa.transdecoder.genome.gff3"
#     output: "results/transdecoder/transcripts.fa.transdecoder.genome.non-red.gff3"
#     run:
#         with open(input[0], 'r') as infile, open(output[0], 'w') as out:
#             prec = ''
#             for line in infile:
#                 if len(line.strip().split('\t'))< 3:
#                     continue
#                 feature_type = line.strip().split('\t')[2]
#                 if feature_type not in ['gene', 'mRNA']:
#                     continue

#                 if prec == feature_type:
#                     continue

#                 out.write(line)

#                 prec = feature_type

rule select_mikado_train:
    input: "results/mikado/mikado.loci.metrics.tsv", "results/mikado/mikado.loci.gff3"
    output: "results/augustus_train/training.gff3"
    conda: "../envs/python3.yaml"
    shell: "python scripts/select_mik_train.py -f 0.5 -e 2 {input} -o {output}"


rule train_augustus:
    input: gff = "results/augustus_train/training.gff3",
           g = config["genome"]+".masked"
    output: "results/augustus_training/autoAugTrain/training/test/augustus.2.withoutCRF.out"
    conda: '../envs/augustus.yaml'
    params: spname = GENOME, odir = "results/augustus_training"
    shell: "mkdir -p {params.odir} && autoAugTrain.pl --trainingset={input.gff} --genome={input.g} "
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
    shell: "python ../scripts/mikado_gtf_to_hints.py {input} -o {output}"


# rule metaeuk_sim_prot_db:

#     "cat *.fasta | mmseqs createdb stdin ${dbname}_DB "
# "mmseqs cluster {dbname}_DB {dbname}_clust_DB tmp/ "
# "mmseqs createsubdb {dbname}_clust_DB {dbname}_DB {dbname}_RepDb "
# "mmseqs createsubdb {dbname}_clust_DB {dbname}_DB_h {dbname}_RepDb_h"
# "mmseqs result2profile {dbname}_RepDb {dbname}_DB {dbname}_clust_DB {dbname}ProfileDb"


rule metauk_genome_db:
    input: config["genome"]+'.masked'
    output: f"results/metaeuk/{GENOME}.DB"
    shell: "metaeuk createdb {input} {output}"

rule metaeuk:
    input: genome =  f"results/metaeuk/{GENOME}.DB", profile = config["metaeuk_profileDB"]
    output: "results/metaeuk/pred_results.gff"
    params: tmp = "results/metaeuk/tmp/"
    threads: 30
    shell: "metaeuk easy-predict {input.genome} {input.profile} {output} {params.tmp} --threads {threads}"

rule metaeuk_to_hint:
    input: "results/metaeuk/pred_results.gff"
    output: "results/hints_for_augustus/protsim_hints.gff"
    conda: "../envs/python3.yaml"
    shell: "python ../scripts/metaeuk_to_hints.py {input} -o {output}"

rule augustus:
    input: h1 = "results/hints_for_augustus/intronhints.gff",
           h2 = "results/hints_for_augustus/mikado_exons_hints.gff",
           h3 = "results/hints_for_augustus/protsim_hints.gff",
           c = config.get("augustus_config", "../../config/extrinsic.E.W.P.cfg"),
           genome = config["genome"] + ".masked"
    output: f"results/augustus/{GENOME}.aug.out"
    params: sp = GENOME
    conda: '../envs/augustus.yaml'
    shell: "augustus --uniqueGeneId=true --gff3=on --species={params.sp} --hintsfile={input.hints} "#to do probably cat hint files
           "--extrinsicCfgFile={input.c} --allow_hinted_splicesites=atac "
           "--alternatives-from-evidence=false {input.genome} > {output}"
