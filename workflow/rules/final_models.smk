"""
Final busco and pfam, pep file and gtf final
"""

if config['metaeuk_only']:
	rule meta_to_gtf:
		input: 
			pasa = "results/filter_models/gene_models.filt.0.6.gff3",
		output: f"results/final_annotation/{GENOME}.gtf"
		params: config["species_3letters_code"] 
		script: "../scripts/format_final_gtf.py"

else:
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
	shell: "transeq {input} {output} && sed -i 's/_1$//g' {output}"


rule get_cds_longest_transcript:
	input: f"results/final_annotation/{GENOME}-cds_all.fa"
	output: f"results/final_annotation/{GENOME}-cds_longest-isoform.fa"
	script: "../scripts/cds_longest_transcript.py"


rule get_pep_longest_transcript:
	input: f"results/final_annotation/{GENOME}-cds_longest-isoform.fa"
	output: f"results/final_annotation/{GENOME}-pep_longest-isoform.fa"
	conda: "../envs/emboss.yaml"
	shell: "transeq {input} {output} && sed -i 's/_1$//g' {output}"


rule blast_swissprot:
	input: db = config["diamond"], sp = f"results/final_annotation/{GENOME}-pep_longest-isoform.fa"
	output: f"results/final_annotation/swissprot_blast_result.blp"
	conda: "../envs/mikado.yaml"
	threads: 5
	shell: "diamond blastp --query {input.sp} --db {input.db} --max-hsps 1 --evalue 1e-5 --outfmt 6 --out {output} -p {threads}"


rule add_genenames:
	input: gtf = f"results/final_annotation/{GENOME}.gtf", blasts = f"results/final_annotation/swissprot_blast_result.blp"
	output: f"results/final_annotation/{GENOME}.with.names.gtf", "results/final_annotation/table_ids_names.tsv"
	script: "../scripts/add_names_swissblast.py"


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


rule make_bed:
	input: f"results/final_annotation/{GENOME}.with.names.gtf"
	output: temp(f"results/final_annotation/{GENOME}.with.names.tmp.bed")
	conda: "../envs/bedops.yaml"
	shell: "gtf2bed < {input} > {output}"


rule list_longest_transcripts:
	input: f"results/final_annotation/{GENOME}-pep_longest-isoform.fa"
	output: f"results/final_annotation/{GENOME}_longest-isoform.txt"
	shell: "grep '>' {input} | sed 's/>//g' > {output}"


rule format_bed:
	input: lgest = f"results/final_annotation/{GENOME}_longest-isoform.txt", bed = f"results/final_annotation/{GENOME}.with.names.tmp.bed"
	output: f"results/final_annotation/{GENOME}.with.names.bed"
	shell: """
	grep -Fwf {input.lgest} {input.bed} | grep mRNA | cut -f 1,2,3,4,6,10 | cut -f 1,6 -d '"' | sed 's/gene_id "//g' | awk -F '\t' '{{$(NF-1)="." FS $(NF-1);}}1' OFS='\t'> {output}
	"""