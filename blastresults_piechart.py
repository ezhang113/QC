#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#handles command line input parameter
import sys

larp6_file = sys.argv[1]

fig, ax = plt.subplots()

df = pd.read_csv(larp6_file, sep='\t')
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

#print(df)
#parsing through tsv file to select only sseqid column


#use groupby to find frequency of each sseq, stored in new dataframe
sseq_count_series = df['sseqid'].value_counts()

#value_counts returns a pandas series so convert to a data frame
sseq_count_df = pd.DataFrame({'sseqid':sseq_count_series.index, 'count':sseq_count_series.values})

#colors = plt.get_cmap('plasma')
count = sseq_count_df['count']
sseq = sseq_count_df['sseqid']

ax.pie(count, labels = sseq, colors=None, radius=3, center=(4, 4),
        wedgeprops={"linewidth": 1, "edgecolor": "white"}, frame=True)


plt.title('sseq counts summary LARP6')
plt.show(block=True)



