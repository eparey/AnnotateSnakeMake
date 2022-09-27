#!/usr/bin/env python

import sys
from pybedtools import BedTool
from collections import defaultdict

def gffKeptEVM(x):
    if ('Parent' in x.attrs and x.attrs['Parent'] in skept) or ('ID' in x.attrs and x.attrs['ID'] in skept) or ('Name' in x.attrs and x.attrs['Name'] in skept):
        return True

input_gtf = snakemake.input[0] #"../GenomeAnnotationAmphiura/results/augustus_train/training.gff3" #
input_bed = snakemake.input[1] #"../GenomeAnnotationAmphiura/resources/genome/Afil_fr2py.fa.out.bed" #
output = snakemake.output[0] #'test.gff' #


cutoff_f=0.1
cutoff_F=0.5

genes = BedTool(input_gtf)
genes = genes.remove_invalid().saveas()

exons=genes.filter(lambda x: x[2] == 'CDS').saveas()
exons=genes.filter(lambda x: x[2] == 'exon').saveas()
print(len(exons), 'exons...')

nbExGene=defaultdict(int)
for i, exon in enumerate(exons):
    nbExGene[exon.attrs['Parent']]+=1

repeats=BedTool(input_bed)

exofilt=exons.intersect(repeats, f=cutoff_f, v=True)

nbFiltEx=defaultdict(int)

for exon in exofilt:
    nbFiltEx[exon.attrs['Parent']]+=1

kept=[]
histex=defaultdict(int)
for gene in nbExGene:
    nfilt=nbFiltEx.get(gene, 0)
    nexon=nbExGene[gene]
    histex['{0}:{1}'.format(nfilt, nexon)]+=1
    if nfilt/float(nexon)>0.9:
        kept.append(gene)

print('initial number of genes:', len(nbExGene))
print('number of genes after filtering', len(kept))

skept=set(kept)


filt_genes=genes.filter(gffKeptEVM).saveas(output)

# with open(output, 'w') as outl:
#     outl.write('\n'.join(list(kept)))