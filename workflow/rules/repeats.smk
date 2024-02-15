
rule make_db:
    """
    Create database for RepeatModeler
    """
    input: config["genome"]
    output: touch('.db_rm_ok')
    log: "logs/repeats/make_db.log"
    params: name = lambda w, input: Path(input[0]).stem
    conda: "../envs/repeats.yaml"
    shell: "BuildDatabase -name {params.name} -engine ncbi {input} 2>&1 | tee {log}"

# rule repeat_mod:
#     """
#     Run RepeatModeler
#     # input: expand(f'{GENOME}.{{ext}}', ext=['nsq', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'translation'])

#     """
#     input: '.db_rm_ok'
#     output: tmp1 = temp(f'{GENOME}-families.fa'), tmp2 = temp(f'{GENOME}-families.stk'),
#             final1 = f"results/repeats/{GENOME}-families.fa", final2 = f"results/repeats/{GENOME}-families.stk"
#     log: "logs/repeats/repeat_modeler.log"
#     params: name = GENOME
#     conda: "../envs/repeats.yaml"
#     threads: min(workflow.cores/4, 12) #RM uses four times the specified number of cores
#     shell: "RepeatModeler -engine ncbi -pa {threads} -database {params.name} 2>&1 | tee {log} && rm -rf RM_* && "
#            "cp {output.tmp1} {output.final1} && cp {output.tmp2} {output.final2}"

rule repeat_masker:
    """
    Run RepeatMasker using the de novo RepeatModeler library
    """
    input: lib = f"results/repeats/{GENOME}-families.fa", g = config["genome"]
    output: config["genome"]+'.masked', config["genome"]+'.align', config["genome"]+'.out', config["genome"]+'.tbl'
    log: "logs/repeats/repeat_masker.log"
    conda: "../envs/repeats.yaml"
    threads: min(workflow.cores/4, 16) #RM uses four times the specified number of cores
    shell: "RepeatMasker -pa {threads} -xsmall -gff -lib {input.lib} {input.g} -a 2>&1 | tee {log}"


rule repeat_masker_to_bed:
    input: config["genome"]+'.out'
    output: f'results/repeats/{GENOME}_repeats.bed'
    shell:"""
          awk -v OFS="\t" 'NR>3,$1=$1' {input} | cut -f 5,6,7,10 > {output}
          """

rule prepare_landscape:
    """
    Prepare inputs to create the repeat landscape
    """
    input: config["genome"]+'.align'
    output: f'results/repeats/{GENOME}_summary.divsum'
    conda: "../envs/repeats.yaml"
    log: "logs/repeats/prepare_landscape.log"
    shell: "calcDivergenceFromAlign.pl -s {output} {input} 2>&1 | tee {log}"


rule repeat_landscape:
    """
    Plot the repeat landscape
    """
    input: div = f'results/repeats/{GENOME}_summary.divsum', tbl = config["genome"]+'.tbl'
    output: f'results/repeats/{GENOME}_repeat_landscape.html'
    log: "logs/repeats/repeat_landscape.log"
    params: gsize = lambda w, input: get_genome_size(input.tbl)
    conda: "../envs/repeats.yaml" 
    shell: "createRepeatLandscape.pl -div {input.div} -g {params.gsize} > {output} 2>&1 | tee {log}"

