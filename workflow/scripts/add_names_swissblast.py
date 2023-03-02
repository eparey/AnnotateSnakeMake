#!/usr/bin/env python

import sys, csv
from collections import defaultdict 

def parseAtt(attributes):
    attrib=dict((att.split()[0], att.split()[1].strip('\"')) for att in attributes.rstrip(';').split(';'))
    return attrib

n=0
gtfFile = snakemake.input[0]
bltFile = snakemake.input[1]
outgtf = snakemake.output[0]
outtable = snakemake.output[1]

blastelts = [(elts[0], elts[1].split('|')[-1]) for elts in csv.reader(open(bltFile), delimiter='\t')]

blast, nnames = {}, defaultdict(int)

for qwr,sub in blastelts:
    #print(qwr,sub)
    gid = qwr.rsplit('.',1)[0]
    fn = sub.split('_')[0]
    fn = fn[0]+fn[1:].lower()
    if not gid in blast:
        #fn=names[sid]
        blast[gid] = fn
        nnames[fn] += 1

gene_names={}
knw_nb=defaultdict(int)
unchr_nb=0
attList=['gene_id','transcript_id','gene_name','gene_type','exon_number']
print(len(blast),'genes with blast hit!')
with open(snakemake.output[0],'w') as out:
    for line in open(gtfFile):
        if not line.strip() or line.startswith('#'): continue
        scaf, source, feat, start, end, score, strand, phase, attributes = line.rstrip().split('\t')
        if feat in ('gene','transcript'): continue
        #print(line)
        attrs=parseAtt(attributes)
        gid=attrs['gene_id']
        #Yprint gid,gene_names.get(gid,'Not found')
        if not gid in gene_names:
            if gid in blast:
                hit=blast[gid]
                ngn=nnames[hit]
                if ngn==1:
                    gene_names[gid]=hit
                elif ngn>1:
                    knw_nb[hit]+=1
                    gene_names[gid]="{0}-{1}".format(hit,knw_nb[hit])
            else:
                unchr_nb+=1
                gene_names[gid]="Unchar_{0}".format(unchr_nb)
        gname=gene_names[gid]
        attrs['gene_name']=gname
        fattrib= ' '.join(["{0} \"{1}\";".format(att, attrs[att]) for att in attList if att in attrs])
        out.write('\t'.join([scaf, source, feat, start, end, score, strand, phase, fattrib])+'\n')

with open(snakemake.output[1],'w') as log: 
    for g in gene_names:    
        #desc=descs.get(blast[g],['','None'])[1] if g in blast else 'None'
        log.write("{}\t{}\n".format(g, gene_names[g]))