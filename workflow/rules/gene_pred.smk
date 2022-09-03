
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
    input: "mikado.loci.metrics.tsv", "mikado.loci.gff3"
    output: "training.gff3"
    shell: "python scripts/select_mik_train.py -f 0.5 -e 2 {input}"

rule train_augustus:
    input: gff = "training.gff3",
           g = config["genome"]+".masked"
    output: directory("results/augustus")
    conda: '../envs/augustus.yaml'
    params: spname = GENOME
    shell: "mkdir -p 'results/augustus' && autoAugTrain.pl --trainingset={input.gff} --genome={input.g} "
           "--species={params.spname} --workingdir={output} --optrounds=1 --verbose --useexisting"


rule portcullis_to_hint:
    input: "results/portcullis/portcullis.filtered.bam"
    output: "results/hints_for_augustus/intronhints.gff"
    conda: '../envs/augustus.yaml'
    shell: "bam2hints --intronsonly --minintronlen=15 --in={input} --out={output}"


rule mikado_to_hint:
    "python gtfToHintsMik.py mikado.loci.gtf"


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
    output: "results/metaeuk/pred_results"
    params: tmp = "results/metaeuk/tmp/"
    threads: 30
    shell: "metaeuk easy-predict {input.genome} {input.profile} {output} {params.tmp} --threads {threads}"

rule metaeuk_to_hint:

rule augustus:
"augustus --uniqueGeneId=true --gff3=on \
   --species=$sp \
   --hintsfile=$hints \
   --extrinsicCfgFile=extrinsic.E.W.P.cfg \
   --allow_hinted_splicesites=atac \
   --alternatives-from-evidence=false $FA > ${FA%.*}.aug.out"
