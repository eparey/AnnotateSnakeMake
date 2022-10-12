import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



# infiles = ["results/pfam/unfiltered.genes.pfam-domains.txt"] + [f"results/pfam/filt{i}.genes.pfam-domains.txt" for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]]
# labels = ["no_filter"] + [f"filter_{i}" for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]]


infiles = list(snakemake.input)
labels = str(snakemake.params).split(",")
output = snakemake.output[0]

rec = []
for i, infile in enumerate(infiles):
    df = pd.read_csv(infile, comment='#', delim_whitespace=True,
                     names=["gene", "a", "b", "c", "domain_id", "domain_name", "type", "e", "f", "g", "h", "i", "j", "k"])
    df = df[["gene", "domain_id", "domain_name", "type"]]
    l = labels[i]
    r = len([i for i in df["type"] if i == "Repeat"])
    d = len(df["domain_id"].unique())

    rec.append([l, r, d])

df = pd.DataFrame.from_records(rec, columns=["gene_set", "total_pfam_hits_repeat", "unique_pfam_hits"])

dfm = df.melt('gene_set', var_name='pfam', value_name='pfam domains')

plot = sns.catplot(x="gene_set", y="pfam domains", hue='pfam', data=dfm, kind='point')

for axes in plot.axes.flat:
    _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=45, ha='right')
plt.tight_layout()
plt.savefig(output)
plt.close('all')