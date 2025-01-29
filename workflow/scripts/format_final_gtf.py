
pasa_gff = snakemake.input[0]
outfile = snakemake.output[0]
sp_code = snakemake.params[0]


def parse_att_gff(attr):
    return dict(at.split('=') for at in attr.strip(';').split(';'))

        
gene_nb, tr_nb, tr_tot = 0, 0, 0 
gID = ''
no_gene_feat = False
with open(pasa_gff, "r") as f, open(outfile, 'w') as out:
    for line in f:
        if not line.strip() or line.startswith('#'): continue
        scaf, source, feat, start, end, score, strand, phase, attributes = line.rstrip().split('\t')
        attrib = parse_att_gff(attributes)

        if feat=='gene':
            gene_nb += 1
            tr_nb = 0
            gID="{0}{1:05d}".format(sp_code, gene_nb)
            attr="gene_id \"{0}\";".format(gID,)
            out.write('\t'.join([scaf, 'AHP', feat, start, end, score, strand, phase, attr])+'\n')

        if feat=='mRNA' or feat=='transcript':


            if not gID or no_gene_feat == True:
                no_gene_feat = True
                gene_nb += 1
                tr_nb = 0
                gID="{0}{1:05d}".format(sp_code, gene_nb)
                attr="gene_id \"{0}\";".format(gID,)
                out.write('\t'.join([scaf, 'AHP', 'gene', start, end, score, strand, phase, attr])+'\n')

            tr_nb += 1
            tr_tot += 1
            trID="{0}{1:05d}.{2}".format(sp_code, gene_nb, tr_nb)
            attr="gene_id \"{0}\"; transcript_id \"{1}\";".format(gID, trID)
            out.write('\t'.join([scaf, 'AHP', feat, start, end, score, strand, phase, attr])+'\n')


        if feat=='exon' or feat=='CDS' or feat=='five_prime_UTR' or feat=='three_prime_UTR':
            attr="gene_id \"{0}\"; transcript_id \"{1}\";".format(gID, trID)
            out.write('\t'.join([scaf, 'AHP', feat, start, end, score, strand, phase, attr])+'\n')


print(gene_nb,'genes, ', tr_tot,'transcripts.')