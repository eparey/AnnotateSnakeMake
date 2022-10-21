
infile = snakemake.input[0]
outfile = snakemake.output[0]

print(infile)

res = {}
seq, gene = '', ''
with open(infile, 'r') as f:

    for line in f:

        if line[0] == '>':

            if seq:
                if gene not in res:
                    res[gene] = (seq, transcript)
                elif len(seq) > len(res[gene][0]):
                    res[gene] = (seq, transcript)
                seq = ''

            transcript = line[1:-1]
            gene = '.'.join(transcript.split('.')[:-1])

        else:
            seq += line


    if seq:
        if gene not in res:
            res[gene] = (seq, transcript)
        elif len(seq) > len(res[gene][0]):
            res[gene] = (seq, transcript)



with open(outfile, 'w') as out:
    for gene in res:
        seq, transcript = res[gene]
        out.write(f'>{transcript}\n')
        out.write(seq)