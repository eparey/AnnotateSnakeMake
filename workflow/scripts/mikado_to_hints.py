#!/usr/bin/env python

import sys
import itertools
import argparse

def parseAtt(attributes):
    attrib=dict(att.split() for att in attributes.rstrip(';').split(';'))
    return attrib


parser = argparse.ArgumentParser()
parser.add_argument("-o", default='mikado_hints.gff', help="output hint file")
args = parser.parse_args()

hinType='mikado'

with open(args.o, 'w') as out, open(sys.argv[1], 'r') as infile:
    for line in infile:
        if line.startswith('#') or line=='\n': continue
        (scaf, source, feat, start, end, score, strand, phase, attributes)=line.rstrip().split('\t')
        if feat == 'exon':
            att = parseAtt(attributes)
            hatt = "group={0};pri=4;source=E".format(att['transcript_id'])
            otl = [scaf.split(';')[0], hinType, 'exonpart', start, end, score, strand, phase, hatt]
            out.write('\t'.join(otl)+'\n')