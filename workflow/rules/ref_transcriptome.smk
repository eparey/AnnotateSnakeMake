
rule star_index:
    """
    Build a genome index for STAR
    """
    input: config["genome"]+".masked"
    output: directory("resource/"+jobname+"/star_index")
    log: "logs/"+jobname+"/star_index_genome.log"
    conda: "../envs/star.yaml"
    threads: 6
    shell: "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} "
           "--genomeFastaFiles {input} 2> {log}"


def get_fq(sample_sheet):
    return


rule star_align:
    """
    Align reads to the genome using STAR
    """
    input:
        fq = get_fq(config["samples"]) #make smart use of wildcards here 
        index = "resource/"+jobname+"/star_index",
    output: "results/"+jobname+"/star/{sample}/Aligned.sortedByCoord.out.bam"
    log: "logs/"+jobname+"/star/{sample}.log"
    params: pref = "results/"+jobname+"/star/" + lambda wildcards: wildcards.sample + "/"
    conda: "../envs/star.yaml"
    threads: 24
    shell: "STAR --genomeDir {input.index} --runThreadN {threads} --readFilesIn {input.fq} "
           "--readFilesCommand zcat --outFileNamePrefix {params.pref} "
           "--outSAMtype BAM SortedByCoordinate --twopassMode Basic 2> {log}"


rule stringtie_assemble:
    """
    Assemble the reference-based transcriptomes with StringTie
    """
    input: "results/"+jobname+"/star/{sample}/Aligned.sortedByCoord.out.bam"
    output: "results/"+jobname+"/stringtie/{sample}_assembly.gtf"
    log: "logs/"+jobname+"/stringtie/{sample}.log"
    params: "--rf"
    conda: "../envs/stringtie.yaml"
    threads: 4
    shell: "stringtie -o {output} {params} -p {threads} {input} 2> {log}"


rule list_stringtie_assemblies:
    """
    Get all stringtie assemblies that we want to merge with taco
    """
    input: infiles = expand("results/"+jobname+"/stringtie/{sample}_assembly.gtf", sample=SAMPLE)
    output: outfile = "results/"+jobname+"/stringtie/assemblies_gtf.txt"
    run:
        with open(output.out, 'w') as out:
            for filename in infiles:
                out.write(filename+'\n')


rule taco_merge:
    """
    Make a consensus transcriptome from the different transcriptome assemblies with taco
    """
    input: "results/"+jobname+"/stringtie/assemblies_gtf.txt"
    output: "results/"+jobname+"/taco/transcripts.gtf"
    log: "logs/"+jobname+"/taco/taco.log"
    params: odir = "results/"+jobname+"/taco/"
    conda: "../envs/taco.yaml"
    threads: 8
    shell: "taco_run -p {threads} -o {params.odir} {input.str_gtf} 2> {log}"


rule transcripts_to_fasta:
    """
    Prepare a fasta file of transcripts for TransDecoder
    """
    input: t = "results/"+jobname+"/taco/transcripts.gtf", g = config["genome"]+".masked"
    output: "results/"+jobname+"/trans_decoder/transcripts.fa"
    log: "logs/"+jobname+"/transdecoder/to_fasta.log"
    conda: "../envs/trans_decoder.yaml"
    shell: "util/gtf_genome_to_cdna_fasta.pl {input.t} {input.g} > {output} 2> {log}"


rule transcripts_to_gff3:
    """
    Prepare a gff3 file of transcripts for TransDecoder
    """
    input: "results/"+jobname+"/taco/transcripts.gtf"
    output: "results/"+jobname+"/trans_decoder/transcripts.gff3"
    log: "logs/"+jobname+"/transdecoder/to_gff3.log"
    conda: "../envs/trans_decoder.yaml"
    shell: "util/gtf_to_alignment_gff3.pl {input} > {output} 2> {log}"


rule trans_decoder:
    """
    Find coding regions within transcripts with TransDecoder
    """
    input:
        fa = "results/"+jobname+"/trans_decoder/transcripts.fa",
        gff = "results/"+jobname+"/trans_decoder/transcripts.gff3"
    output:
        gff = "results/"+jobname+"/trans_decoder/transcripts.fa.transdecoder.genome.gff3",
        bed = "results/"+jobname+"/trans_decoder/transcripts.fa.transdecoder.genome.bed"
    log:
    conda: "../envs/trans_decoder.yaml"
    shell:
        "TransDecoder.LongOrfs -t {input.tfa} && "
        "util/cdna_alignment_orf_to_genome_orf.pl {input.fa}.transdecoder.gff3 {input.gff} {input.fa} > {output.gff} && "
        "utilgff3_file_to_bed.pl {output.gff} > {output.bed}"
