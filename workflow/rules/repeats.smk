from pathlib import Path

GENOME = Path(config["genome"]).stem

rule make_db:
    """
    Create database for RepeatModeler
    """
    input: config["genome"]
    output: temp(expand(f'{GENOME}.{{ext}}', ext=['nsq', 'nhr', 'nin', 'nnd', 'nni', 'nog', 'translation']))
    log: "logs/repeats/make_db.log"
    params: name = lambda w, input: Path(input).stem
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
    output: config["genome"]+'.masked', config["genome"]+'.align'
    log: "logs/repeats/repeat_masker.log"
    conda: "../envs/repeats.yaml"
    threads: 12 #will use three or four times this I think
    shell: "RepeatMasker -pa {threads} -xsmall -gff -lib {input.lib} {input.g} -a 2>&1 | tee {log}"


#will need a post-activate script to fix repeatmasker I think ?
#To use `calcDivergenceFromAlign.pl` from RepeatMasker with a conda installation:
#--> need to open the script and modify the hard-coded path to their perl modules... Replace 'use lib "/home/rhubley/projects/RepeatMasker";' line 88 by 'use lib "/import/kg_dyows01/workspace1/parey/miniconda3/envs/repeat_mod/share/RepeatMasker";' <-- should be $CONDA_PREFIX/share/RepeatMasker


#TODO: make repeatlandscape as part of the generated report??
#rule prepare_landscape
#rule repeat_landscape