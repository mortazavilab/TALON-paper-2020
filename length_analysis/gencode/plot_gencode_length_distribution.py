from Bio import SeqIO
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

infile = "../../refs/GENCODE_v29/gencode.v29.transcripts.fa"


dists = []
for record in SeqIO.parse(infile, "fasta"):
    dists.append(len(record))

# Plot histogram
fname = "gencode_v29_isoform_length_distribution.png"
x = pd.Series(dists, name="Isoform length (kb)")
x = x/1000.0
plt.xlim(0,20)
ax = sns.distplot(x, kde = False, color = "dodgerblue", 
                  bins = np.arange(0, max(x), 1))#1000)

# Plot vertical lines at interesting points
style = dict(size=10, color='black')
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

med = round(np.median(x), 1)
quantile_90 = round(x.quantile(.90), 1)
quantile_99 = round(x.quantile(.99), 1)

plt.axvline(med, linestyle = '--', color = 'lightgrey')
plt.axvline(quantile_90, linestyle = '--', color = 'lightgrey')
plt.axvline(quantile_99, linestyle = '--', color = 'lightgrey')

ax.text(med + 0.1, ymax*7/8, "50%: " + str(med) + " kb", **style)
ax.text(quantile_90 + 0.1, ymax*6.5/8, "90%: " + str(quantile_90) + " kb", **style)
ax.text(quantile_99 + 0.1, ymax*6/8, "99%: " + str(quantile_99) + " kb", **style)

plt.title("Length distribution of GENCODE v29 isoforms")
plt.tight_layout()
plt.savefig(fname, dpi = 600, bbox_inches='tight')
plt.close()
