rule trinity:
    input: fq1 = [SAMPLES[i][0] for i in SAMPLES],
           fq2 = [SAMPLES[i][1] for i in SAMPLES],
           samples = config["rna_samples"]
    output: "results/trinity/Trinity.fasta"
    threads: 24
    params: odir = lambda w, output: os.path.dirname(output[0]), cpu_bfly = 4
    # conda: "../envs/trinity.yaml"
    shell: "Trinity --trimmomatic --seqType fq --samples_file {input.samples} --CPU {threads} "
           "--max_memory 300G --output {params.odir}  --full_cleanup"



rule gmap_db:
    input: config["genome"] + ".masked"
    output: directory("results/gmapdb")
    params: g = GENOME
    conda: "../envs/gmap.yaml"
    shell: "gmap_build -D {output} -d {params.g} {input}"


rule gmap:
    input: fa = "results/trinity/Trinity.fasta", db = directory("results/gmapdb/")
    output: 'results/trinity/Trinity_transcripts.gff3'
    params: g = GENOME
    threads: 20
    shell: "gmap -d {params.g} -D {input.db} -f 3 -n 0 -x 50 -t {threads} -B 5 --gff3-add-separators=0 "
           "--intronlength=500000 {input.fa} > {output}"


#minimap2 -t 16 -x splice:hq -uf -c $GNM $TRS > PriCau_trin.paf
#paftools.js  splice2bed PriCau_trin.paf >  PriCau_trin.bed


# Trinity --CPU 24 --max_memory 500G \
#     --seqType fq --full_cleanup \
#     --trimmomatic \
#     --left $LRD --right $RRD --output $SFX\_trinity