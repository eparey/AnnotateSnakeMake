#!/usr/bin/env python

genelist_file = snakemake.input[1] 
input_gff = snakemake.input[0]
output = snakemake.output[0]


with open(genelist_file, 'r') as infile:
    genelist = {line.strip() for line in infile}

with open(input_gff, 'r') as infile, open(output, 'w') as out:
    for line in infile:
        if line[0] == '#' or not line.strip():
            continue
        att = line.strip().split('\t')[8]
        feat = line.strip().split('\t')[2]
        ID = att[3:].split(';')[0]

        # if feat == "gene" and ID in genelist:
        #     out.write(line)
        if '.'.join(ID.split('.')[:2]) in genelist:
            out.write(line)
