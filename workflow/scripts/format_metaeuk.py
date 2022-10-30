#!/usr/bin/env python

input_gtf = snakemake.input[0]
output = snakemake.output[0]


with open(input_gtf, 'r') as infile, open(output, 'w') as out:
    for line in infile:
        line_split = line.strip().split('\t')
        feat = line_split[2]
        attr = {entry.split('=')[0]:entry.split('=')[1] for entry in line_split[-1].split(';')}

        if feat == 'gene':
            attr['ID'] = attr['Target_ID']
            attr.pop('Target_ID')
            attr.pop('TCS_ID')

        else:
            attr['ID'] = attr['Target_ID'] + '_' + feat + attr['TCS_ID'].split(feat)[-1]
            attr['Parent'] = attr['Target_ID']
            if feat in ['CDS', 'exon']:
                attr['Parent'] += '_mRNA'
            attr.pop('Target_ID')
            attr.pop('TCS_ID')


        attr_write = ''
        for key in attr:
            attr_write += key + '=' + attr[key] + ';'

        line_split[-1] = attr_write

        out.write('\t'.join(line_split)+'\n')