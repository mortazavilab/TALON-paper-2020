import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('sirv_table.tsv', sep='\t')
df.columns = ['Technology', 'Known', 'Novel', 'ISM']
df['Novel minus ISMs'] = df['Novel'] - df['ISM']

df = pd.melt(df, id_vars='Technology', var_name='isoform_category',
             value_vars=['Known', 'Novel', 'Novel minus ISMs'],
             value_name='num_isoforms')
df.rename({'num_isoforms': 'Number of detected isoforms'},
           axis=1, inplace=True)
colors = {'Known': "#009E73", 'Novel': "#56B4E9", 'Novel minus ISMs': "#E69F00"}
ax = sns.barplot(x='Technology',
                 y='Number of detected isoforms',
                 hue='isoform_category',
                 data=df,
                 palette=colors,
                saturation=1
                )
ax.set(xlim=(-0.5,3.5))
ax.plot([-0.5, 3.5], [69,69], color='gray', linestyle='--')
plt.legend(bbox_to_anchor=(1.05, 0.6), loc=2, borderaxespad=0.)
font = {'color': 'gray',
       'size': 12}
plt.text(3.55, 67.5, "69 known SIRV isoforms", fontdict=font)
plt.savefig('sirv_detection_by_technology.png', bbox_inches='tight')
