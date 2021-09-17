PE_ORDER = """Initial bases num
Initial reads num
Read 1 Total written (filtered) Round 1
Read 1 Trimmed bases Round 1
Read 1 basepairs processed Round 1
Read 1 with adapter Round 1
Read 1 with adapter percent Round 1
Read 2 Total written (filtered) Round 1
Read 2 Trimmed bases Round 1
Read 2 basepairs processed Round 1
Read 2 with adapter Round 1
Read 2 with adapter percent Round 1
Reads after cutadapt 1
Reads Written percent Round 1
Reads that were too short percent Round 1
Too short reads Round 1
Total written (filtered) Round 1
Total written (filtered) percent Round 1
Trimmed bases Round 1
Trimmed bases percent Round 1
Processed bases Round 2
Processed reads Round 2
Read 1 Total written (filtered) Round 2
Read 1 Trimmed bases Round 2
Read 1 basepairs processed Round 2
Read 1 with adapter Round 2
Read 1 with adapter percent Round 2
Read 2 Total written (filtered) Round 2
Read 2 Trimmed bases Round 2
Read 2 basepairs processed Round 2
Read 2 with adapter Round 2
Read 2 with adapter percent Round 2
Reads after cutadapt 2
Reads Written percent Round 2
Reads that were too short percent Round 2
Too short reads Round 2
Total written (filtered) Round 2
Total written (filtered) percent Round 2
Trimmed bases Round 2
Trimmed bases percent Round 2
Percent Repetitive
Repetitive Reads
STAR genome input reads
% of reads mapped to multiple loci
% of reads mapped to too many loci
% of reads unmapped: other
% of reads unmapped: too many mismatches
% of reads unmapped: too short
Average input read length
Average mapped length
Deletion average length
Deletion rate per base
Insertion average length
Insertion rate per base
Mismatch rate per base, percent
Number of reads mapped to multiple loci
Number of reads mapped to too many loci
Number of splices: AT/AC
Number of splices: Annotated (sjdb)
Number of splices: GC/AG
Number of splices: GT/AG
Number of splices: Non-canonical
Number of splices: Total
STAR genome uniquely mapped %
STAR genome uniquely mapped
Usable reads
removed_count
total_count
Clipper peaks num
Percent usable / mapped
Percent Usable / Input""".split("\n")

SE_ORDER = """Initial bases num
Initial reads num
Reads with adapter Round 1
Reads with adapter percent Round 1
Reads after cutadapt 1
Reads Written percent Round 1
Reads that were too short percent Round 1
Too short reads Round 1
Total written (filtered) Round 1
Total written (filtered) percent Round 1
Trimmed bases Round 1
Trimmed bases percent Round 1
Processed bases Round 2
Processed reads Round 2
Reads with adapter Round 2
Reads with adapter percent Round 2
Reads after cutadapt 2
Reads Written percent Round 2
Reads that were too short percent Round 2
Too short reads Round 2
Total written (filtered) Round 2
Total written (filtered) percent Round 2
Trimmed bases Round 2
Trimmed bases percent Round 2
Percent Repetitive
Repetitive Reads
STAR genome input reads
% of reads mapped to multiple loci
% of reads mapped to too many loci
% of reads unmapped: other
% of reads unmapped: too many mismatches
% of reads unmapped: too short
Average input read length
Average mapped length
Deletion average length
Deletion rate per base
Insertion average length
Insertion rate per base
Mismatch rate per base, percent
Number of reads mapped to multiple loci
Number of reads mapped to too many loci
Number of splices: AT/AC
Number of splices: Annotated (sjdb)
Number of splices: GC/AG
Number of splices: GT/AG
Number of splices: Non-canonical
Number of splices: Total
STAR genome uniquely mapped %
STAR genome uniquely mapped
Usable reads
Clipper peaks num
Percent usable / mapped
Percent Usable / Input""".split("\n")



slim_qc_metrics = ["Initial reads num", "Reads after cutadapt 2", "Repetitive Reads", "Percent Repetitive", "STAR genome input reads",
                   "STAR genome uniquely mapped", "STAR genome uniquely mapped %",
                   'Number of reads mapped to too many loci',
                   '% of reads unmapped: too short', '% of reads mapped to too many loci', "Usable reads",
                   "Percent usable / mapped", "Percent Usable / Input", "Clipper peaks num"]

