#select intepreter and change environment to anaconda
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

plt.style.use('_mpl-gallery-nogrid')
fig, ax = plt.subplots()

#parsing through tsv file
#ideally use a loop to loop through all the files in that directory --> figure out how to connect this with bash
for larp6_file in ():
    df = pd.read_csv(larp6_file, sep='\t', usecols = ['sseqid'])
    #use groupby to find frequency of each sseq, stored in new dataframe 
    sseq_count = df.groupby('sseqid').count()

    colors = plt.get_cmap(plt.Colormap)

    ax.pie(sseq_count['count'], labels = sseq_count['sseqid'], colors=colors, radius=3, center=(4, 4),
       wedgeprops={"linewidth": 1, "edgecolor": "white"}, frame=True)
    ax.set()

    plt.title('sseq counts...')
    plt.show()


#df = pd.read_csv(larp6_file, sep='\t', usecols = ['sseqid'])

#load only sseqid results from TSV files into a data frame
#df = pd.read_csv('LARP6_CTRL_IN2_rand10p.unmappedblast', sep='\t', usecols = ['sseqid'])


