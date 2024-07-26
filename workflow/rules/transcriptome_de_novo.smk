rule trinity:
    input: fq1 = [SAMPLES[i][0] for i in SAMPLES],
           fq2 = [SAMPLES[i][1] for i in SAMPLES],
           samples = config["rna_samples"]
    output: "results/trinity/Trinity.fasta"
    threads: 24
    params: odir = lambda w, output: os.path.dirname(output[0]) + '/'
    conda: "../envs/trinity.yaml"
    log: "logs/trinity/trinity.log"
    shell: "Trinity --trimmomatic --seqType fq --samples_file {input.samples} --CPU {threads} "
           "--max_memory 300G --output {params.odir} --full_cleanup 2>&1 | tee {log} && "
           "mv results/trinity.Trinity.fasta {output} && mv results/trinity.Trinity.fasta.gene_trans_map {params.odir}"

rule gmap_db:
    input: config["genome"] + ".masked"
    output: directory("results/gmapdb")
    params: g = GENOME
    conda: "../envs/gmap.yaml"
    log: "logs/trinity/gmap_db.log"
    shell: "gmap_build -D {output} -d {params.g} {input} 2>&1 | tee {log}"


rule gmap:
    input: fa = "results/trinity/Trinity.fasta", db = "results/gmapdb/"
    output: 'results/trinity/Trinity_transcripts.gff3'
    params: g = GENOME
    conda: "../envs/gmap.yaml"
    threads: 20
    log: "logs/trinity/gmap.log"
    shell: "gmap -d {params.g} -D {input.db} -f 3 -n 0 -x 50 -t {threads} -B 5 --gff3-add-separators=0 "
           "--intronlength=500000 {input.fa} > {output}"


#minimap2 -t 16 -x splice:hq -uf -c $GNM $TRS > trin.paf
#paftools.js  splice2bed trin.paf >  trin.bed