#!/usr/bin/env python

import sys
import itertools
import argparse

def parseAtt(attributes):
    attrib=dict(att.split("=") for att in attributes.rstrip(';').split(';'))
    return attrib


# parser = argparse.ArgumentParser()
# parser.add_argument("-i", help="input gtf file")
# parser.add_argument("-o", default='mikado_hints.gff', help="output hint file")
# parser.add_argument("-t", default='mikado', help="hint type")
# args = parser.parse_args()

hinType=snakemake.params[0]
source=snakemake.params[1]
featType=snakemake.params[2]

with open(snakemake.output[0], 'w') as out, open(snakemake.input[0], 'r') as infile:
    for line in infile:
        if line.startswith('#') or line=='\n': continue
        (scaf, source, feat, start, end, score, strand, phase, attributes)=line.rstrip().split('\t')
        if feat == 'exon':
            att = parseAtt(attributes)
            hatt = "group={0};pri=4;source={1}".format(att['Parent'], source)
            otl = [scaf.split(';')[0], hinType, featType, start, end, score, strand, phase, hatt]
            out.write('\t'.join(otl)+'\n')