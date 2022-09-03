#!/usr/bin/env python

import sys
import itertools

def parseAtt(attributes):
    attrib=dict(att.split() for att in attributes.rstrip(';').split(';'))
    return attrib

hinType='mikado'

out=open(sys.argv[1].rsplit('.', 1)[0]+'.exh.gff','w')
for line in open(sys.argv[1]):
    if line.startswith('#') or line=='\n': continue
    (scaf,source,feat,start,end,score,strand,phase,attributes)=line.rstrip().split('\t')
    if feat=='exon':
        att=parseAtt(attributes)
        hatt="group={0};pri=4;source=E".format(att['transcript_id'])
        otl=[scaf.split(';')[0],hinType,'exonpart',start,end,score,strand,phase,hatt]
        out.write('\t'.join(otl)+'\n')