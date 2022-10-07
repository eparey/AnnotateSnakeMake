
aug_file = snakemake.input[0]
aug_file_ok = snakemake.output[0]

with open(aug_file, 'r') as infile, open(aug_file_ok, 'w') as out:
    for line in infile:
        feature = line.strip().split('\t')[2]
        if feature == 'CDS':
            out.write(line)
            exon_line = line.replace("CDS", "exon").replace("cds", "exon")
            out.write(exon_line)

        if feature == 'gene':
            out.write(line)

        if feature == 'transcript':
            rna_line = line.replace("transcript", "mRNA")
            out.write(rna_line)


