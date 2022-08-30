
rule star_index:
    """
    Build a genome index for STAR
    """
    input: config["genome"]+".masked"
    output: directory("resources/star_index")
    log: "logs/star_index_genome.log"
    conda: "../envs/star.yaml"
    threads: 6
    shell: "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} "
           "--genomeFastaFiles {input} 2>&1 | tee {log}"


# def get_fq(sample_sheet):
#     return
SAMPLE = ["1d-M-1_P6959_2024_subset"]

rule star_align:
    """
    Align reads to the genome using STAR
    """
    input: fq1 = "resources/{sample}_R1_tr.fq.gz", fq2 = "resources/{sample}_R2_tr.fq.gz", index = "resources/star_index"
    output: "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    log: "logs/star/{sample}.log"
    params: pref =  lambda wildcards, input: wildcards.sample
    conda: "../envs/star.yaml"
    threads: 24
    shell: "STAR --genomeDir {input.index} --runThreadN {threads} --readFilesIn {input.fq1} {input.fq2} "
           "--readFilesCommand zcat --outFileNamePrefix results/star/{params.pref}/ "
           "--outSAMtype BAM SortedByCoordinate --twopassMode Basic 2>&1 | tee {log}"


rule stringtie_assemble:
    """
    Assemble the reference-based transcriptomes with StringTie
    """
    input: "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output: "results/stringtie/{sample}_assembly.gtf"
    log: "logs/stringtie/{sample}.log"
    params: "--rf"
    conda: "../envs/stringtie.yaml"
    threads: 4
    shell: "stringtie -o {output} {params} -p {threads} {input} 2>&1 | tee {log}"


rule list_stringtie_assemblies:
    """
    Get all stringtie assemblies that we want to merge with taco
    """
    input: infiles = expand("results/stringtie/{sample}_assembly.gtf", sample=SAMPLE)
    output: outfile = "results/stringtie/assemblies_gtf.txt"
    run:
        with open(output.outfile, 'w') as out:
            for filename in input.infiles:
                out.write(filename+'\n')

rule taco_merge:
    """
    Make a consensus transcriptome from the different transcriptome assemblies with taco
    """
    input: "results/stringtie/assemblies_gtf.txt"
    output: "results/taco/assembly.gtf"
    log: "logs/taco/taco.log"
    params: odir = "results/taco/"
    conda: "../envs/taco.yaml"
    threads: 8
    shell: "rm -r {params.odir} && taco_run -p {threads} -o {params.odir} {input} 2>&1 | tee {log}"


rule transcripts_to_fasta:
    """
    Prepare a fasta file of transcripts for TransDecoder
    """
    input: t = "results/taco/assembly.gtf", g = config["genome"]+".masked"
    output: "results/trans_decoder/transcripts.fa"
    log: "logs/transdecoder/to_fasta.log"
    conda: "../envs/trans_decoder.yaml"
    shell: "gtf_genome_to_cdna_fasta.pl {input.t} {input.g} > {output} 2>&1 | tee {log}"


rule transcripts_to_gff3:
    """
    Prepare a gff3 file of transcripts for TransDecoder
    """
    input: "results/taco/assembly.gtf"
    output: "results/trans_decoder/transcripts.gff3"
    log: "logs/transdecoder/to_gff3.log"
    conda: "../envs/trans_decoder.yaml"
    shell: "gtf_to_alignment_gff3.pl {input} > {output} 2>&1 | tee {log}"


rule transdecoder_predict:
    """
    Find coding regions within transcripts with TransDecoder
    """
    input:
        fa = "results/trans_decoder/transcripts.fa",
        gff = "results/trans_decoder/transcripts.gff3"
    output:
        tmp_gff = "transcripts.fa.transdecoder.gff3"
    log: log1 = "logs/transdecoder/transdecoder_longorfs.log", log2 = "logs/transdecoder/transdecoder_predict.log"
    conda: "../envs/trans_decoder.yaml"
    params: config.get("transdecoder_predict_args", "")
    shell:
        "TransDecoder.LongOrfs -t {input.fa} 2>&1 | tee {log.log1} && "
        "TransDecoder.Predict {params} -t {input.fa} 2>&1 | tee {log.log2}"


rule transdecoder_to_genome:
    input: fa = "results/trans_decoder/transcripts.fa", gff =  "results/trans_decoder/transcripts.gff3", gff_v2 = "transcripts.fa.transdecoder.gff3"
    output: gff = "results/trans_decoder/transcripts.fa.transdecoder.genome.gff3",
            bed = "results/trans_decoder/transcripts.fa.transdecoder.genome.bed"
    conda: "../envs/trans_decoder.yaml"
    log: log1 = "logs/transdecoder/transdecoder.log", log2 = "logs/transdecoder/to_bed.log"
    shell:
        "cdna_alignment_orf_to_genome_orf.pl {input.gff_v2} {input.gff} {input.fa} > {output.gff} 2>&1 | tee {log.log1} && "
        "gff3_file_to_bed.pl {output.gff} > {output.bed} 2>&1 | tee {log.log2} && rm pipeliner* && rm -r transcripts.fa*"