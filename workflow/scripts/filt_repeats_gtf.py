#!/usr/bin/env python

import sys
from pybedtools import BedTool
from collections import defaultdict

def gffKeptEVM(x):
    if x.attrs['gene_id'] in skept:
        return True

input_gtf = snakemake.input[0]
input_bed = snakemake.input[1]
output = snakemake.output[0]


cutoff_f=0.5
cutoff_F=0.5

genes = BedTool(sys.argv[1])
genes = genes.remove_invalid().saveas()

exons=genes.filter(lambda x: x[2] == 'CDS').saveas()
#exons=genes.filter(lambda x: x[2] == 'exon').saveas()
print(len(exons), 'exons...')

nbExGene=defaultdict(int)
for i, exon in enumerate(exons):
    nbExGene[exon.attrs['gene_id']]+=1

repeats=BedTool(sys.argv[2])

exofilt=exons.intersect(repeats, f=cutoff_f, v=True)

nbFiltEx=defaultdict(int)

for exon in exofilt:
    nbFiltEx[exon.attrs['gene_id']]+=1

kept=[]
histex=defaultdict(int)
for gene in nbExGene:
    nfilt=nbFiltEx.get(gene, 0)
    nexon=nbExGene[gene]
    histex['{0}:{1}'.format(nfilt, nexon)]+=1
    if nfilt/float(nexon)>0.5:
        kept.append(gene)

print('initial number of genes:', len(nbExGene))
print('number of genes after filtering', len(kept))

skept=set(kept)


filt_genes=genes.filter(gffKeptEVM).saveas(outf)

with open(output, 'w') as outl:
    outl.write('\n'.join(list(kept)))