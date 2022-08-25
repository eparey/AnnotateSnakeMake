#will need a post-activate script to fix repeatmasker I think ?
#To use `calcDivergenceFromAlign.pl` from RepeatMasker with a conda installation:
#--> need to open the script and modify the hard-coded path to their perl modules... Replace 'use lib "/home/rhubley/projects/RepeatMasker";' line 88 by 'use lib "/import/kg_dyows01/workspace1/parey/miniconda3/envs/repeat_mod/share/RepeatMasker";' <-- should be $CONDA_PREFIX/share/RepeatMasker

#will need to check RM output files path

GENOME = os.path.basename(config["genome"])

rule make_db:
    """
    """
    input: config["genome"]
    output:
    log: "logs/"+jobname+"/repeats/make_db.log"
    params: name = GENOME
    conda: "../envs/repeats.yaml"
    shell: "BuildDatabase -name {params.name} -engine ncbi {input}"

rule repeat_mod:
    input: 
    output: GENOME+"-families.fa"
    log: "logs/"+jobname+"/repeats/repeat_modeler.log"
    params: name = GENOME
    conda: "../envs/repeats.yaml"
    threads: 12 #will use three or four times this I think
    shell: "RepeatModeler -engine ncbi -pa {threads} -database {params.name}"

rule repeat_masker:
    input: GENOME+"-families.fa"
    output: 
    log: "logs/"+jobname+"/repeats/repeat_masker.log"
    params: name = GENOME
    conda: "../envs/repeats.yaml"
    threads: 12 #will use three or four times this I think
    shell: "RepeatMasker -pa {threads} -xsmall -gff -lib {input.lib} {input.genome} -a {output.ali}"

#TODO: make repeatlandscape as part of the generated report??
#rule prepare_landscape
#rule repeat_landscape

