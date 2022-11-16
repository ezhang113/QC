#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
import xmltodict
from xml.dom import minidom


#handles command line input parameter
import sys

larp6_file = sys.argv[1]

fig, ax = plt.subplots()

df = pd.read_csv(larp6_file, sep='\t')
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

#use groupby to find frequency of each sseq, stored in new dataframe
sseq_count_series = df['sseqid'].value_counts()

#loop through series to determine which elements to remove and add into "other" column
to_remove = []
other_count = 0

for index,values in sseq_count_series.iteritems():
        if(values < 50000):
                to_remove.append(index)
                other_count += values

print(sseq_count_series.size)

#remove these elements from series
sseq_count_series2 = sseq_count_series.drop(to_remove)
print(sseq_count_series2.size)

#replace index sseqids with ncbi name
for index in sseq_count_series2.iteritems():
        sseq_count_series.rename(index={index:#call ncbi query function})

print(sseq_count_series)

#generate new pandas series with new element to concatenate with old series
#d = {'Other':other_count}
#ser = pd.Series(data=d; index=['Other']

#append new element
#sseq_count = sseq_count_series2.append(ser)

#value_counts returns a pandas series so convert to a data frame
sseq_count_df = pd.DataFrame({'sseqid':sseq_count_series.index, 'count':sseq_count_series.values})

#count = sseq_count_df['count']
#sseq = sseq_count_df['sseqid']

#ax.pie(count, labels = sseq, colors=None, radius=3, center=(4, 4),
#	wedgeprops={"linewidth": 1, "edgecolor": "white"}, frame=True)


#plt.title('sseq counts summary LARP6')
#plt.show(block=True)
#plt.savefig("sseq_counts1")

