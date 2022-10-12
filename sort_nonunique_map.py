#extract non-uniquely-mapped reads from *genome-mapped.bam file & convert into 
# a format that BLAST understand

#Pysam: python module that allows you to read and manipulate data stored in 
#SAM/BAM files

#to figure out the multipmapped reads
import pysam #figure out how to import this

#read file in BAM format, AlignmentFile object
CTRL_IP1_genomemapped = pysam.AlignmentFile('LARP6.CTRL_IN1.umi.r1.fq.genome-mapped.bam', 'rb')
#b for compressed BAM u for uncompressed, figure out which one

#fetch() method reads single sequence and returns AlignedSegment object
for sequence in CTRL_IP1_genomemapped.fetch(){
    print(read)
    #where you determine where multimapped reads are aligning to
}

#to align sequences to BLAST by converting to FASTA format 



