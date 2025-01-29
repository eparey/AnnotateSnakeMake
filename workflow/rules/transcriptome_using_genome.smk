
sample_file = config.get("rna_samples", "")
# print(sample_file)
SAMPLES = {}
if sample_file:
    with open(sample_file, 'r') as infile:
        for line in infile:
            if line[0] != '#':
                sample, rep, fastq1, fastq2 = line.strip().split('\t')
                assert fastq1.endswith('.gz') and fastq2.endswith('.gz'), "Error: expected gzipped fastqs, STAR command will likely not run properly."
                SAMPLES[sample] = [fastq1, fastq2]

# print(SAMPLES)

rule star_index:
    """
    Build a genome index for STAR
    """
    input: f"{GENOME_PATH}.masked"
    output: temp(directory("resources/star_index"))
    log: "logs/star_index_genome.log"
    conda: "../envs/star.yaml"
    threads: 6
    shell: "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} "
           "--genomeFastaFiles {input} 2>&1 | tee {log}"


rule star_align:
    """
    Align reads to the genome using STAR
    """
    input: fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
           fq2 = lambda wildcards: SAMPLES[wildcards.sample][1],
           index = "resources/star_index"
    output: "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    log: "logs/star/{sample}.log"
    params: pref =  lambda wildcards, input: wildcards.sample
    conda: "../envs/star.yaml"
    threads: 24
    shell: "STAR --genomeDir {input.index} --runThreadN {threads} --readFilesIn {input.fq1} {input.fq2} "
           "--readFilesCommand zcat --outFileNamePrefix results/star/{params.pref}/ --outSAMstrandField intronMotif "
           "--outSAMtype BAM SortedByCoordinate --twopassMode Basic 2>&1 | tee {log}"


rule stringtie_assemble:
    """
    Assemble the reference-based transcriptomes with StringTie
    """
    input: "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output: "results/stringtie/{sample}_assembly.gtf"
    log: "logs/stringtie/{sample}.log"
    conda: "../envs/stringtie.yaml"
    threads: 4
    shell: "stringtie -o {output} -p {threads} {input} 2>&1 | tee {log}"


rule list_stringtie_assemblies:
    """
    Get all stringtie assemblies that we want to merge with taco
    """
    input: infiles = expand("results/stringtie/{sample}_assembly.gtf", sample=list(SAMPLES.keys()))
    output: outfile = "results/stringtie/assemblies_gtf.txt"
    log: "logs/taco/prepare.log"
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
    params: odir = lambda w, output: os.path.dirname(output[0])
    conda: "../envs/taco.yaml"
    threads: 8
    shell: "rm -r {params.odir} && taco_run -p {threads} -o {params.odir} {input} 2>&1 | tee {log}"


rule merge_bams:
    """
    Merge bams for portcullis
    """
    input: expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output: "results/star/merged.bams"
    conda: '../envs/samtools.yaml'
    shell: "samtools merge {output} {input}"


rule portcullis:
    """
    Extract curated splice junctions
    """
    input: b = "results/star/merged.bams", g = config["genome"]+'.masked'
    output: "results/portcullis/3-filt/portcullis_filtered.pass.junctions.bed", temp("results/portcullis/portcullis.filtered.bam")
    params: odir = "results/portcullis/"
    conda: '../envs/portcullis.yaml'
    threads: 20
    shell: "portcullis full -t {threads} -v --bam_filter -o {params.odir} {input.g} {input.b}"
