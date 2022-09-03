#!/usr/bin/env python

import sys, csv
import argparse
from collections import defaultdict

def parseAtt(attr):
    return dict(at.split('=') for at in attr.rstrip(';').split(';'))

parser = argparse.ArgumentParser()
parser.add_argument("metrics", help="mikado metrics file")
parser.add_argument("loci", help="mikado loci gff3 file")
parser.add_argument("-o", default='training.gff3', help="output training gff3 file")
parser.add_argument("-f", default=0.5, type=float, help="fraction of cds covered by blast (default:0.5)")
parser.add_argument("-e", default=2, type=int, help="minimum number of exons (default:2)")
args = parser.parse_args()


#Usage: select_mik_train.py mikado.loci.metrics.tsv mikado.loci.gff3
if not len(sys.argv) >= 3:
    sys.exit('Usage: select_mik_train.py <mikado.loci.metrics.tsv> <mikado.loci.gff3> ')

training=set()
genes=set()
with open(args.metrics) as mtr:
    for loc in csv.reader(mtr, delimiter='\t'):
        if loc[0]=='tid': continue
        tid, gid, blastCov, hasStart, hasStop, nbExon = loc[0], loc[2], float(loc[7]), loc[29], loc[30], int(loc[25])
        if blastCov>args.f and hasStart=='True' and hasStop=='True' and nbExon>=args.e:
            if not gid in genes:
                training.add(tid)
                genes.add(gid)
            #print gid,blastCov,hasStart,hasStop,nbExon

print "{0} transcrits and {1} genes selected!".format(len(training),len(genes))
print list(training)[0:5]

added=set()
with open(args.o, 'w') as tgf:
    for line in open(args.loci):
        if not line.strip() or line.startswith('#'): continue
        scaf, source, feat, start, end, score, strand, phase, attributes = line.rstrip().split('\t')
        att=parseAtt(attributes)
        if feat=='mRNA':
            print att['ID'],att['ID'] in training
            if att['ID'] in training:
                tgf.write(line)
        elif feat=='exon' or feat=='CDS':
            if att['Parent'] in training:
                tgf.write(line)
