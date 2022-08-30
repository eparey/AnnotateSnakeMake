from pathlib import Path

GENOME = Path(config["genome"]).stem

rule make_db:
    """
    Create database for RepeatModeler
    """
    input: config["genome"]
    output: temp(expand(f'{GENOME}.{{ext}}', ext=['nsq', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'translation']))
    log: "logs/repeats/make_db.log"
    params: name = lambda w, input: Path(input[0]).stem
    conda: "../envs/repeats.yaml"
    shell: "BuildDatabase -name {params.name} -engine ncbi {input} 2>&1 | tee {log}"

rule repeat_mod:
    """
    Run RepeatModeler
    """
    input: expand(f'{GENOME}.{{ext}}', ext=['nsq', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'translation'])
    output: f'{GENOME}-families.fa', f'{GENOME}-families.stk'
    log: "logs/repeats/repeat_modeler.log"
    params: name = lambda w, input: Path(input[0]).stem
    conda: "../envs/repeats.yaml"
    threads: 12 #will use three or four times this I think
    shell: "RepeatModeler -engine ncbi -pa {threads} -database {params.name} 2>&1 | tee {log} && rm -rf RM_*"

rule mv_consensus:
    input: fa = f'{GENOME}-families.fa', stk = f'{GENOME}-families.stk'
    output: fa = f"resources/{GENOME}-families.fa", stk = f"resources/{GENOME}-families.stk"
    log: "logs/repeats/mv_consensus.log"
    shell: "mv {input.fa} {output.fa} && mv {input.stk} {output.stk} 2>&1 | tee {log}"

rule repeat_masker:
    """
    Run RepeatMasker using the de novo RepeatModeler library
    """
    input: lib = f"resources/{GENOME}-families.fa", g = config["genome"]
    output: config["genome"]+'.masked', config["genome"]+'.align', config["genome"]+'.tbl'
    log: "logs/repeats/repeat_masker.log"
    conda: "../envs/repeats.yaml"
    threads: 12 #will use three or four times this I think
    shell: "RepeatMasker -pa {threads} -xsmall -gff -lib {input.lib} {input.g} -a 2>&1 | tee {log}"

rule prepare_landscape:
    input: config["genome"]+'.align'
    output: config["genome"]+'_summary.divsum'
    conda: "../envs/repeats.yaml"
    shell: "calcDivergenceFromAlign.pl -s {output} {input}"

def get_genome_size(input_file):
    with open(input_file, 'r') as infile:
        for line in infile:
            if "total length:" in line:
                size = int(line.split()[2])
                return size

rule repeat_landscape:
    input: div = config["genome"]+'_summary.divsum', tbl = config["genome"]+'.tbl'
    output: config["genome"]+'_repeat_landscape.html'
    params: gsize = lambda w, input: get_genome_size(input.tbl)
    conda: "../envs/repeats.yaml" 
    shell: "createRepeatLandscape.pl -div {input.div} -g {params.gsize} > {output}"

#TODO: add plot to report