rule mikado_configure:
    input:
        junctions = "results/portcullis/3-filt/portcullis_filtered.pass.junctions.bed",
        fa = config["diamond"].replace(".dmnd", ".fa"),
        g = config["genome"]+'.masked',
        assemblies = config["mikado"]
    output: "results/mikado/configuration.yaml", "results/mikado/human_scoring.yaml"
    params: odir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/mikado.yaml"
    shell:"mikado configure --list {input.assemblies} --reference {input.g}  --mode permissive "
          "--scoring human.yaml --copy-scoring {output[1]} --junctions {input.junctions} "
          "-bt {input.fa} --out-dir {params.odir} {output[0]}"


rule mikado_prepare:
    input: "results/mikado/configuration.yaml"
    output: "results/mikado/mikado_prepared.fasta",  "results/mikado/mikado_prepared.gtf"
    conda: "../envs/mikado.yaml"
    shell: "mikado prepare --json-conf {input}"


rule mikado_diamond:
    input: q = "mikado_prepared.fasta", db = config["diamond"]
    output: "results/mikado/mikado_diamond.xml"
    threads: 20
    conda: "../envs/mikado.yaml"
    shell: "diamond blastx --query {input.q} --max-target-seqs 5 --sensitive --index-chunks 1 "
           "--threads {threads} --db {input.db} --evalue 1e-6 --outfmt 5 --out {output}"


rule transcripts_to_gff3:
    """
    Prepare a gff3 file of transcripts for TransDecoder
    """
    input: "results/mikado/mikado_prepared.gtf"
    output: "results/mikado/mikado_prepared.gff3"
    log: "logs/transdecoder/to_gff3.log"
    conda: "../envs/transdecoder.yaml"
    shell: "gtf_to_alignment_gff3.pl {input} > {output} 2>&1 | tee {log}"


rule transdecoder_predict:
    """
    Find coding regions within transcripts with TransDecoder
    """
    input:
        fa = "results/mikado/mikado_prepared.fasta",
        gff = "results/mikado/mikado_prepared.gff3"
    output:
        tmp_gff = "transcripts.fa.transdecoder.gff3"
    log: log1 = "logs/transdecoder/transdecoder_longorfs.log", log2 = "logs/transdecoder/transdecoder_predict.log"
    conda: "../envs/transdecoder.yaml"
    params: config.get("transdecoder_predict_args", "")
    shell:
        "TransDecoder.LongOrfs -t {input.fa} 2>&1 | tee {log.log1} && "
        "TransDecoder.Predict {params} -t {input.fa} 2>&1 | tee {log.log2}"


rule transdecoder_to_genome:
    input: fa = "results/mikado/mikado_prepared.fasta", gff =  "results/mikado/mikado_prepared.gff3", gff_v2 = "transcripts.fa.transdecoder.gff3"
    output: gff = "results/mikado/mikado_transdecoder.gff3",
            bed = "results/mikado/mikado_transdecoder.bed"
    conda: "../envs/transdecoder.yaml"
    log: log1 = "logs/transdecoder/transdecoder.log", log2 = "logs/transdecoder/to_bed.log"
    shell:
        "cdna_alignment_orf_to_genome_orf.pl {input.gff_v2} {input.gff} {input.fa} > {output.gff} 2>&1 | tee {log.log1} && "
        "gff3_file_to_bed.pl {output.gff} > {output.bed} 2>&1 | tee {log.log2} && rm pipeliner* && rm -r transcripts.fa*"


rule mikado_serialize:
    input: conf = "results/mikado/configuration.yaml", blast = "results/mikado/mikado_diamond.xml",
           fa = config["diamond"].replace(".dmnd", ".fa"), bed = "results/mikado/mikado_transdecoder.bed"
    output: "results/mikado/mikado.subloci.gff3"
    conda: "../envs/mikado.yaml"
    shell: "mikado serialise --json-conf {input.conf} --xml {input.blast} --orfs {input.bed} "
           "--blast_targets {input.fa}"

rule mikado_pick:
    input: "results/mikado/mikado.subloci.gff3"
    output: "results/mikado/mikado.loci.metrics.tsv", "results/mikado/mikado.loci.gff3"
    conda: "../envs/mikado.yaml"
    shell: "mikado pick --json-conf configuration.yaml --subloci_out {input}"


rule mikado_to_gtf:
    input: "results/mikado/mikado.loci.gff3"
    output: "results/mikado/mikado.loci.gtf"
    conda: "../envs/agat.yaml"
    shell: "agat_convert_sp_gxf2gxf.pl --gff {input}"