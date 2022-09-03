rule trinity:
    input: fq1 = [SAMPLES[i][0] for i in SAMPLES],
           fq2 = [SAMPLES[i][1] for i in SAMPLES]
    output: "results/trinity/Trinity.fasta"
    threads: 70
    params: odir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/trinity.yaml"
    shell: "Trinity --seqType fq --left {input.fq1} --right {input.fq2} --CPU {threads} --max_memory 100G --output {params.odir}"


rule gmap_db:
    input: config["genome"] + ".masked"
    output: directory("results/gmapdb")
    params: g = GENOME
    conda: "../envs/gmap.yaml"
    shell: "gmap_build -D {output} -d {params.g} {input}"


rule gmap:
    input: fa = "results/trinity/Trinity.fasta", db = directory("results/gmapdb/")
    output: 'results/trinity_mapped/Trinity_transcripts.gff3'
    params: g = GENOME
    threads: 20
    shell: "gmap -d {params.g} -D {input.db} -f 3 -n 0 -x 50 -t {threads} -B 5 --gff3-add-separators=0 "
           "--intronlength=500000 {input.fa} > {output}"