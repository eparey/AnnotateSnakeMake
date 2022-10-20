"""
Final busco and pfam, pep file and gtf final
"""


rule to_gtf:
	input: 
		pasa = "results/pasa/pasa_gene_models.gff3",
		clm = "results/pasa/.end"
	output: f"results/final_annotation/{GENOME}.gtf"
	params: config["species_3letters_code"] 
	script: "../scripts/format_final_gtf.py"


rule get_cds:
	input: 
		gtf = f"results/final_annotation/{GENOME}.gtf",
		genome = config["genome"] + ".masked"
	output: f"results/final_annotation/{GENOME}-cds_all.fa"
	conda: "../envs/pasa.yaml"
	shell: "gffread -x {output} -g {input.genome} {input.gtf}"


rule get_pep:
	input: f"results/final_annotation/{GENOME}-cds_all.fa"
	output: f"results/final_annotation/{GENOME}-pep_all.fa"
	conda: "../envs/emboss.yaml"
	shell: "transeq {input} {output} && sed -i 's/_1$//' {output}"


rule get_cds_longest_transcript:
	input: f"results/final_annotation/{GENOME}-cds_all.fa"
	output: f"results/final_annotation/{GENOME}-cds_longest-isoform.fa"
	script: "../scripts/cds_longest_transcript.py"


rule get_pep_longest_transcript:
	input: f"results/final_annotation/{GENOME}-cds_longest-isoform.fa"
	output: f"results/final_annotation/{GENOME}-pep_longest-isoform.fa"
	conda: "../envs/emboss.yaml"
	shell: "transeq {input} {output} && sed -i 's/_1$//' {output}"


rule busco_final:
    input: f"results/final_annotation/{GENOME}-pep_longest-isoform.fa"
    output: "results/final_annotation/busco/short_summary.specific.metazoa_odb10.busco.json"
    params: jname = "busco", odir = "results/final_annotation/"
    threads: 4
    conda: "../envs/busco.yaml"
    shell: "busco -l metazoa_odb10 --tar --mode proteins -o {params.jname} -f -i {input} --cpu {threads} --out_path {params.odir}"


rule pfam_final:
    input: f = f"results/final_annotation/{GENOME}-pep_longest-isoform.fa", d = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/Pfam-A.hmm.h3i"
    output: "results/final_annotation/pfam/pfam-domains.txt"
    conda: "../envs/pfam.yaml"
    params: d = "/home/elise/projects/annot/AnnotateSnakeMake/resources/pfam_db/"
    threads: 8
    shell: "pfam_scan.pl -fasta {input.f} -dir {params.d} -outfile {output} -cpu {threads}"