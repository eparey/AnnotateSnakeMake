import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

import json

infiles = list(snakemake.input)
labels = str(snakemake.params).split(",")
output = snakemake.output[0]

buscos = {l:f for l, f in zip(labels, infiles)}

cols = ["Complete", "Multi copy", "Single copy", "Missing", "Fragmented"]
cols += ["Complete percentage", "Multi copy percentage", "Single copy percentage", "Missing percentage", "Fragmented percentage"]

d_res = {}
for busco in buscos:
    infile = buscos[busco]
    with open(infile, 'r') as j:
        res = json.loads(j.read())
    d_res[busco] = {k:v for k, v in res["results"].items() if k in cols}


df = pd.DataFrame.from_dict(d_res, orient='index')
df['gene_set'] = df.index
dfm = df.melt('gene_set', var_name='busco metazoa', value_name='busco %')
plot = sns.catplot(x="gene_set", y="busco %", hue='busco metazoa', data=dfm, kind='point')

for axes in plot.axes.flat:
    _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=45, ha='right')
plt.tight_layout()
plt.savefig(output)
plt.close('all')