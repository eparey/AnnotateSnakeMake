

rule mikado_configure:
    input:
        junctions = "results/portcullis/3-filt/portcullis_filtered.pass.junctions.bed",
        fa = config["diamond"],
        g = config["genome"]+'.masked',
        assemblies = config["mikado"],
        trinity = 'results/trinity/Trinity_transcripts.gff3',
        taco = "results/taco/assembly.gtf"
    output: temp("results/mikado/configuration.yaml"), "results/mikado/human_scoring.yaml"
    params: odir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/mikado.yaml"
    shell:"mikado configure --list {input.assemblies} --reference {input.g}  --mode permissive "
          "--scoring human.yaml --copy-scoring {output[1]} --junctions {input.junctions} "
          "-bt {input.fa} -od {params.odir} {output[0]}"

rule fix_config:
    input: "results/mikado/configuration.yaml"
    output: "results/mikado/configuration.ok.yaml"
    shell: "sed 's#mikado_prepared.fasta#results/mikado/mikado_prepared.fasta#g' {input} | "
           "sed 's#mikado_prepared.gtf#results/mikado/mikado_prepared.gtf#g' > {output} | "
           "sed 's#/mikado.db#results/mikado/mikado.db#g' > {output}"


rule mikado_prepare:
    input: "results/mikado/configuration.ok.yaml"
    output: "results/mikado/mikado_prepared.fasta",  "results/mikado/mikado_prepared.gtf"
    params: odir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/mikado.yaml"
    shell: "mikado prepare --json-conf {input} --out {output[1]} --out_fasta {output[0]}"

rule diamond_db:
    input:  config["diamond"]
    output:  config["diamond"].replace(".fa", ".dmnd")
    threads: 30
    conda: "../envs/mikado.yaml"
    shell: "diamond makedb --in {input} -d {output}"

rule mikado_diamond:
    input: q = "results/mikado/mikado_prepared.fasta", db = config["diamond"].replace(".fa", ".dmnd")
    output: "results/mikado/mikado_diamond.xml"
    threads: 30
    conda: "../envs/mikado.yaml"
    shell: "diamond blastx --query {input.q} --max-target-seqs 5 --sensitive --index-chunks 1 "
           "--threads {threads} --db {input.db} --evalue 1e-6 --outfmt 5 --out {output}"


rule transdecoder_predict:
    """
    Find coding regions within transcripts with TransDecoder
    """
    input:
        fa = "results/mikado/mikado_prepared.fasta"
    output:
        bed = "results/mikado/mikado_prepared.fasta.transdecoder.bed"
    log: log1 = "logs/transdecoder/transdecoder_longorfs.log", log2 = "logs/transdecoder/transdecoder_predict.log"
    conda: "../envs/transdecoder.yaml"
    params: config.get("transdecoder_predict_args", "")
    shell:
        "TransDecoder.LongOrfs -t {input.fa} 2>&1 | tee {log.log1} && "
        "TransDecoder.Predict {params} -t {input.fa} 2>&1 | tee {log.log2} && "
        "mv mikado_prepared.fasta.transdecoder.bed {output} && "
        "rm pipeliner* && rm -r mikado_prepared.fasta.transdecoder*"


rule mikado_serialize:
    input: conf = "results/mikado/configuration.ok.yaml", blast = "results/mikado/mikado_diamond.xml",
           fa = config["diamond"], bed =  "results/mikado/mikado_prepared.fasta.transdecoder.bed" #"results/mikado/mikado_transdecoder.bed"
    output: "results/mikado/mikado.db"
    params: odir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/mikado.yaml"
    shell: "mikado serialise --json-conf {input.conf} --xml {input.blast} --orfs {input.bed} "
           "--blast_targets {input.fa} -od {params.odir}"


rule mikado_pick:
    input:  c = "results/mikado/configuration.ok.yaml", db="results/mikado/mikado.db" #"mikado.db"
    output: "results/mikado/mikado.loci.metrics.tsv", "results/mikado/mikado.loci.gff3", "results/mikado/mikado.subloci.gff3"
    conda: "../envs/mikado.yaml"
    params: odir = lambda w, output: os.path.dirname(output[0]), subloc = "mikado.subloci.gff3"
    threads: 24
    shell: "mikado pick --json-conf {input.c} --start-method=spawn --subloci_out {params.subloc} -od {params.odir} -db {input.db} --procs={threads}"


rule mikado_to_gtf:
    input: "results/mikado/mikado.loci.gff3"
    output: "results/mikado/mikado.loci.gtf"
    conda: "../envs/agat.yaml"
    shell: "agat_convert_sp_gxf2gxf.pl --gff {input} -o {output}"