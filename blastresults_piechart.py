#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

#handles command line input parameter
import sys

larp6_file = sys.argv[1]

plt.style.use('_mpl-gallery-nogrid')
fig, ax = plt.subplots()

#parsing through tsv file to select only sseqid column
df = pd.read_csv(larp6_file, sep='\t', usecols = ['sseqid'])

#use groupby to find frequency of each sseq, stored in new dataframe 
sseq_count = df.groupby('sseqid').count()

colors = plt.get_cmap(plt.Colormap)

ax.pie(sseq_count['count'], labels = sseq_count['sseqid'], colors=colors, radius=3, center=(4, 4),
    wedgeprops={"linewidth": 1, "edgecolor": "white"}, frame=True)
ax.set()

plt.title('sseq counts summary')
plt.show()




