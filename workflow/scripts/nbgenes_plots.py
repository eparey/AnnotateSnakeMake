import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import json

infiles = list(snakemake.input)
labels = str(snakemake.params).split(",")
output = snakemake.output[0]


for (i, j) in zip(labels, infiles):
    print(i, j)

buscos = {l:f for (l, f) in zip(labels, infiles)}

nb_genes = []
for infile in infiles:
    nb_genes.append(sum(1 for _ in open(infile, 'r')))


rec = [[nb_genes[i], labels[i]] for i in range(len(nb_genes))]
df = pd.DataFrame.from_records(rec, columns=["nb_genes", "gene_set"])

g = sns.barplot(data=df, y="nb_genes", x="gene_set")
sns.despine()
g.set_xticklabels(g.get_xticklabels(), rotation=45, ha='right')

for p in g.patches:
    height = p.get_height()
    g.text(p.get_x()+p.get_width()/2., height + 300, int(height) , ha="center")

plt.tight_layout()

plt.savefig(output)
plt.show()