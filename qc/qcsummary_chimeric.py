#!/usr/bin/env python

"""
General parsers for QC output of pipeline file,
generally pass a handle to the file you want to parse,
returns a dict containing all useful information
Currently this isn't standard
"""

# transitionning to python2/python3 support
# uncomment from this compatibility import list, as py3/py2 support progresses
from __future__ import print_function
from __future__ import division
# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
# from future.utils import iteritems

import glob
import os
import argparse

import pandas as pd
import pybedtools
import pysam

# from parse_cutadapt import parse_cutadapt_file
from .qcsummary_rnaseq import parse_star_file
from qc import parse_cutadapt as pc
from qc import column_names

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from subprocess import call

from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'DejaVu Sans'})

CLASH_SE_ORDER = """Initial bases num
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
STAR repeat input (Chimeric)
STAR repeat uniquely mapped (Chimeric)
STAR repeat mapped to multiple loci (Chimeric)
STAR genome input (Chimeric)
STAR repeat input (Non-chimeric)
STAR repeat uniquely mapped (Non-chimeric)
STAR repeat mapped to multiple loci (Non-chimeric)
STAR genome input (Non-chimeric)
Repeat (Non-chimeric) % of reads mapped to multiple loci
Repeat (Non-chimeric) % of reads mapped to too many loci
Repeat (Non-chimeric) % of reads unmapped: other
Repeat (Non-chimeric) % of reads unmapped: too many mismatches
Repeat (Non-chimeric) % of reads unmapped: too short
Genome (Non-chimeric) % of reads mapped to multiple loci
Genome (Non-chimeric) % of reads mapped to too many loci
Genome (Non-chimeric) % of reads unmapped: other
Genome (Non-chimeric) % of reads unmapped: too many mismatches
Genome (Non-chimeric) % of reads unmapped: too short
Genome (Non-chimeric) Average input read length
Genome (Non-chimeric) Average mapped length
Genome (Non-chimeric) Deletion average length
Genome (Non-chimeric) Deletion rate per base
Genome (Non-chimeric) Insertion average length
Genome (Non-chimeric) Insertion rate per base
Genome (Non-chimeric) Mismatch rate per base, percent
Genome (Non-chimeric) Number of reads mapped to multiple loci
Genome (Non-chimeric) Number of reads mapped to too many loci
Genome (Non-chimeric) Number of splices: AT/AC
Genome (Non-chimeric) Number of splices: Annotated (sjdb)
Genome (Non-chimeric) Number of splices: GC/AG
Genome (Non-chimeric) Number of splices: GT/AG
Genome (Non-chimeric) Number of splices: Non-canonical
Genome (Non-chimeric) Number of splices: Total
Repeat (Chimeric) % of reads mapped to multiple loci
Repeat (Chimeric) % of reads mapped to too many loci
Repeat (Chimeric) % of reads unmapped: other
Repeat (Chimeric) % of reads unmapped: too many mismatches
Repeat (Chimeric) % of reads unmapped: too short
Genome (Chimeric) % of reads mapped to multiple loci
Genome (Chimeric) % of reads mapped to too many loci
Genome (Chimeric) % of reads unmapped: other
Genome (Chimeric) % of reads unmapped: too many mismatches
Genome (Chimeric) % of reads unmapped: too short
Genome (Chimeric) Average input read length
Genome (Chimeric) Average mapped length
Genome (Chimeric) Deletion average length
Genome (Chimeric) Deletion rate per base
Genome (Chimeric) Insertion average length
Genome (Chimeric) Insertion rate per base
Genome (Chimeric) Mismatch rate per base, percent
Genome (Chimeric) Number of reads mapped to multiple loci
Genome (Chimeric) Number of reads mapped to too many loci
Genome (Chimeric) Number of splices: AT/AC
Genome (Chimeric) Number of splices: Annotated (sjdb)
Genome (Chimeric) Number of splices: GC/AG
Genome (Chimeric) Number of splices: GT/AG
Genome (Chimeric) Number of splices: Non-canonical
Genome (Chimeric) Number of splices: Total
STAR genome uniquely mapped (Chimeric)
STAR genome uniquely mapped % (Chimeric)
STAR genome uniquely mapped (Non-chimeric)
STAR genome uniquely mapped % (Non-chimeric)
Usable reads (Chimeric)
Percent Usable / Mapped (Chimeric)
Percent Usable / Input (Chimeric)
Usable reads (Non-chimeric)
Percent Usable / Mapped (Non-chimeric)
Percent Usable / Input (Non-chimeric)
Clipper peaks num (merged)""".split("\n")

clash_slim_qc_metrics = [
    "Initial reads num", 
    "Reads after cutadapt 1", 
    "STAR repeat input (Chimeric)", 
    "STAR repeat uniquely mapped (Chimeric)",
    "STAR repeat mapped to multiple loci (Chimeric)", 
    "STAR genome input (Chimeric)",
    "STAR genome uniquely mapped (Chimeric)", 
    "STAR genome uniquely mapped % (Chimeric)",
    "STAR repeat input (Non-chimeric)", 
    "STAR repeat uniquely mapped (Non-chimeric)",
    "STAR repeat mapped to multiple loci (Non-chimeric)", 
    'Genome (Chimeric) Number of reads mapped to too many loci',
    'Genome (Chimeric) % of reads unmapped: too short', 
    'Genome (Chimeric) % of reads mapped to too many loci', 
    "STAR genome input (Non-chimeric)",
    "STAR genome uniquely mapped (Non-chimeric)", 
    "STAR genome uniquely mapped % (Non-chimeric)",
    'Genome (Non-chimeric) Number of reads mapped to too many loci',
    'Genome (Non-chimeric) % of reads unmapped: too short', 
    'Genome (Non-chimeric) % of reads mapped to too many loci', 
    "Usable reads (Chimeric)",
    "Percent Usable / Mapped (Chimeric)", 
    "Usable reads (Non-chimeric)",
    "Percent Usable / Mapped (Non-chimeric)", 
    "Clipper peaks num (merged)"
]

def write_clipseq_metrics_csv(df, output_csv):
    """
    Writes CSV file with clipseq metrics.
    Writes "annotated" (abridged with the most useful QC stats) file.
    """
    df.to_csv(output_csv)
    df[clash_slim_qc_metrics].to_csv(os.path.splitext(output_csv)[0] + ".annotated.csv")
    
    
def clashseq_metrics(analysis_dir, output_csv, percent_usable, number_usable, cutadapt_round1_suffix, repeat_mapped_suffix, repeat_mapped_eclip_suffix, genome_mapped_suffix, genome_mapped_eclip_suffix, rm_dup_suffix, rm_dup_eclip_suffix, peak_suffix):
    
    df = clipseq_metrics_df(
        analysis_dir=analysis_dir,
        percent_usable=percent_usable,
        number_usable=number_usable,
        cutadapt_round1_suffix=cutadapt_round1_suffix, 
        repeat_mapped_suffix=repeat_mapped_suffix, 
        repeat_mapped_eclip_suffix=repeat_mapped_eclip_suffix,
        genome_mapped_suffix=genome_mapped_suffix, 
        genome_mapped_eclip_suffix=genome_mapped_eclip_suffix,
        rm_dup_suffix=rm_dup_suffix,
        rm_dup_eclip_suffix=rm_dup_eclip_suffix,
        peak_suffix=peak_suffix
    )
    for key in CLASH_SE_ORDER:
        if key not in df.columns:
            print("{} NOT IN DF".format(key))
    df = df[CLASH_SE_ORDER]
    
    write_clipseq_metrics_csv(df, output_csv)
    

def clipseq_metrics_df(
        analysis_dir, percent_usable,
        number_usable,
        cutadapt_round1_suffix, 
        repeat_mapped_suffix, 
        repeat_mapped_eclip_suffix,
        genome_mapped_suffix, 
        genome_mapped_eclip_suffix,
        rm_dup_suffix,
        rm_dup_eclip_suffix,
        peak_suffix,
        num_seps=3,
        sep=".",
    ):
    #######################################
    """
    Reports all clip-seq metrics in a given analysis directory
    outputs must follow gabes naming clipseq pipeline / naming conventions"
    Args:
        analysis_dir:
        num_seps:
        sep:
        percent_usable:
        number_usable:
    Returns:
    """

    cutadapt_round1_names, \
    repeat_mapped_names, \
    repeat_mapped_eclip_names, \
    genome_mapped_names, \
    genome_mapped_eclip_names, \
    rm_duped_names, \
    rm_duped_eclip_names, \
    peaks_names = get_all_names(
        analysis_dir=analysis_dir,
        cutadapt_round1_suffix=cutadapt_round1_suffix,
        repeat_mapped_suffix=repeat_mapped_suffix,
        repeat_mapped_eclip_suffix=repeat_mapped_eclip_suffix,
        genome_mapped_suffix=genome_mapped_suffix, 
        genome_mapped_eclip_suffix=genome_mapped_eclip_suffix,
        rm_dup_suffix=rm_dup_suffix,
        rm_dup_eclip_suffix=rm_dup_eclip_suffix,
        peak_suffix=peak_suffix,
        sep=sep,
        num_seps=num_seps
    )
    ###########################################################################
    # make dataframes for each column
    ###########################################################################
    if len(cutadapt_round1_names) > 0:
        combined_df = pd.DataFrame(
            {
                name: pc.parse_cutadapt_file(cutadapt_file, paired_end=False)
                for name, cutadapt_file in cutadapt_round1_names.items()
            }
        ).transpose()

        combined_df.columns = [
            "{} Round 1".format(col) for col in combined_df.columns
        ]
    if len(repeat_mapped_names) > 0:
        star_rep_df = pd.DataFrame(
            {
                name: parse_star_file(repeat_mapped_file) 
                for name, repeat_mapped_file in repeat_mapped_names.items()
            }
        ).transpose()
        star_rep_df.columns = ['Repeat (Chimeric) {}'.format(x) for x in star_rep_df.columns]
    if len(repeat_mapped_eclip_names) > 0:
        star_rep_eclip_df = pd.DataFrame(
            {
                name: parse_star_file(repeat_mapped_eclip_file) 
                for name, repeat_mapped_eclip_file in repeat_mapped_eclip_names.items()
            }
        ).transpose()
        star_rep_eclip_df.columns = ['Repeat (Non-chimeric) {}'.format(x) for x in star_rep_eclip_df.columns]
    if len(genome_mapped_names) > 0:
        star_genome_df = pd.DataFrame(
            {
                name: parse_star_file(genome_mapped_file) 
                for name, genome_mapped_file in genome_mapped_names.items()
            }
        ).transpose()
        star_genome_df.columns = ['Genome (Chimeric) {}'.format(x) for x in star_genome_df.columns]
    if len(genome_mapped_eclip_names) > 0:
        star_genome_eclip_df = pd.DataFrame(
            {
                name: parse_star_file(genome_mapped_eclip_file) 
                for name, genome_mapped_eclip_file in genome_mapped_eclip_names.items()
            }
        ).transpose()
        star_genome_eclip_df.columns = ['Genome (Non-chimeric) {}'.format(x) for x in star_genome_eclip_df.columns]
    if len(rm_duped_names) > 0:
        rm_duped_df = pd.DataFrame(
            {
                name: parse_rm_duped_metrics_file_se(rm_duped_file, label="Usable reads (Chimeric)")
                for name, rm_duped_file in rm_duped_names.items()
            }
        ).transpose()
    else:
        rm_duped_df = pd.DataFrame(index=['Usable reads (Chimeric)', 'removed_count (Chimeric)', 'total_count (Chimeric)']).T
    if len(rm_duped_eclip_names) > 0:
        rm_duped_eclip_df = pd.DataFrame(
            {
                name: parse_rm_duped_metrics_file_se(rm_duped_eclip_file, label="Usable reads (Non-chimeric)")
                for name, rm_duped_eclip_file in rm_duped_eclip_names.items()
            }
        ).transpose()
    else:
        rm_duped_df = pd.DataFrame(index=['Usable reads (Non-chimeric)', 'removed_count (Non-chimeric)', 'total_count (Non-chimeric)']).T
    if len(peaks_names) > 0:
        peaks_df = pd.DataFrame(
            {
                name: {"Clipper peaks num (merged)": len(pybedtools.BedTool(peaks_file))}
                for name, peaks_file in peaks_names.items()
            }
        ).transpose()
    else:
        peaks_df = pd.DataFrame(index=['Clipper peaks num (merged)']).T
    ###########################################################################

    ###########################################################################
    # get the base dataframe rnaseq metrics dataframe
    ##############################
    ###########################################################################

    ###########################################################################
    # merge dataframes
    ##################
    
    combined_df = pd.merge(combined_df, star_rep_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, star_rep_eclip_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, star_genome_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, star_genome_eclip_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rm_duped_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rm_duped_eclip_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, peaks_df,
                           left_index=True, right_index=True, how="outer")
    
    ###########################################################################
    # Rename columns to be useful
    # combined_df = combined_df.rename(
    #     columns={"Reads Written Round 2": "Reads after cutadapt 2"
    #              })

    ###########################################################################
    print(combined_df.columns)
    combined_df = combined_df.rename(
        columns={
            "Processed bases Round 1": "Initial bases num",
            "Processed reads Round 1": "Initial reads num",
            "Reads Written Round 1": "Reads after cutadapt 1",
            "Repeat (Chimeric) Number of input reads": "STAR repeat input (Chimeric)",
            "Repeat (Chimeric) STAR genome uniquely mapped number": "STAR repeat uniquely mapped (Chimeric)",
            "Repeat (Chimeric) Number of reads mapped to multiple loci": "STAR repeat mapped to multiple loci (Chimeric)",
            "Genome (Chimeric) Number of input reads": "STAR genome input (Chimeric)",
            "Genome (Chimeric) STAR genome uniquely mapped number": "STAR genome uniquely mapped (Chimeric)",
            "Genome (Chimeric) STAR genome uniquely mapped %": "STAR genome uniquely mapped % (Chimeric)",
            
            "Repeat (Non-chimeric) Number of input reads": "STAR repeat input (Non-chimeric)",
            "Repeat (Non-chimeric) STAR genome uniquely mapped number": "STAR repeat uniquely mapped (Non-chimeric)",
            "Repeat (Non-chimeric) Number of reads mapped to multiple loci": "STAR repeat mapped to multiple loci (Non-chimeric)",
            "Genome (Non-chimeric) Number of input reads": "STAR genome input (Non-chimeric)",
            "Genome (Non-chimeric) STAR genome uniquely mapped number": "STAR genome uniquely mapped (Non-chimeric)",
            "Genome (Non-chimeric) STAR genome uniquely mapped %": "STAR genome uniquely mapped % (Non-chimeric)",        
        }
    )
    # compute useful stats
    ######################
    combined_df['STAR genome uniquely mapped (Chimeric)'] = combined_df['STAR genome uniquely mapped (Chimeric)'].astype(float)
    combined_df['Initial reads num'] = combined_df['Initial reads num'].astype(float)
    try:
        combined_df["Percent Usable / Mapped (Chimeric)"] = (combined_df['Usable reads (Chimeric)'] / combined_df['STAR genome uniquely mapped (Chimeric)'])
        combined_df["Percent Usable / Input (Chimeric)"] = (combined_df['Usable reads (Chimeric)'] / combined_df['Initial reads num'])
        combined_df["Percent Usable / Mapped (Non-chimeric)"] = (combined_df['Usable reads (Non-chimeric)'] / combined_df['STAR genome uniquely mapped (Chimeric)'])
        combined_df["Percent Usable / Input (Non-chimeric)"] = (combined_df['Usable reads (Non-chimeric)'] / combined_df['Initial reads num'])
        combined_df["Percent Repetitive"] = 1 - (combined_df['STAR genome input (Chimeric)'] / combined_df['STAR repeat input (Chimeric)'].astype(float))
        combined_df["Repetitive Reads"] = combined_df['STAR repeat input (Chimeric)'] - combined_df['STAR genome input (Chimeric)']
        
    except ZeroDivisionError:
        print("passing on ZeroDivisionError")
        pass

    return combined_df


def get_all_names(
        analysis_dir,
        cutadapt_round1_suffix,
        repeat_mapped_suffix,
        repeat_mapped_eclip_suffix,
        genome_mapped_suffix, 
        genome_mapped_eclip_suffix,
        rm_dup_suffix,
        rm_dup_eclip_suffix,
        peak_suffix,
        sep,
        num_seps
    ):
    ###########################################################################
    # get file paths
    ################

    cutadapt_round1_files = glob.glob(os.path.join(analysis_dir, cutadapt_round1_suffix))
    repeat_mapped_files = glob.glob(os.path.join(analysis_dir, repeat_mapped_suffix))
    repeat_mapped_eclip_files = glob.glob(os.path.join(analysis_dir, repeat_mapped_eclip_suffix))
    genome_mapped_files = glob.glob(os.path.join(analysis_dir, genome_mapped_suffix))
    genome_mapped_eclip_files = glob.glob(os.path.join(analysis_dir, genome_mapped_eclip_suffix))
    rm_duped_files = glob.glob(os.path.join(analysis_dir, rm_dup_suffix))
    rm_duped_eclip_files = glob.glob(os.path.join(analysis_dir, rm_dup_eclip_suffix))
    
    if peak_suffix is not None:
        peaks_files = glob.glob(os.path.join(analysis_dir, peak_suffix))
    ###########################################################################

    ###########################################################################
    # get file names
    ################
    cutadapt_round1_names = get_names(cutadapt_round1_files, num_seps, sep)
    repeat_mapped_names = get_names(repeat_mapped_files, num_seps, sep)
    repeat_mapped_eclip_names = get_names(repeat_mapped_eclip_files, num_seps, sep)
    genome_mapped_names = get_names(genome_mapped_files, num_seps, sep)
    genome_mapped_eclip_names = get_names(genome_mapped_eclip_files, num_seps, sep)
    rm_duped_names = get_names(rm_duped_files, num_seps, sep)
    rm_duped_eclip_names = get_names(rm_duped_eclip_files, num_seps, sep)
    
    if peak_suffix is not None:
        peaks_names = get_names(peaks_files, num_seps, sep)
    else:
        peaks_names = []
    ###########################################################################
    return cutadapt_round1_names, repeat_mapped_names, repeat_mapped_eclip_names, genome_mapped_names, genome_mapped_eclip_names, rm_duped_names, rm_duped_eclip_names, peaks_names


def get_names(files, num_seps, sep):
    """
    Given a list of files,
     return that files base name and the path to that file

    :param files: list
        list of files
    :param num_seps: int
        number of separators to call real names
    :param sep: str
        separator to split names on
    :return basenames: dict
        dict basename to file
    """

    dict_basename_to_file = {
        sep.join(os.path.basename(file).split(sep)[0: num_seps]): file
        for file in files
    }
    # print("get names dict: {}".format(dict_basename_to_file))
    return dict_basename_to_file


def parse_peak_metrics(fn):
    """
    Unused function that has parsed/will parse CLIPPER metrics.

    :param fn: basestring
    :return spot_dict: dict
    """
    with open(fn) as file_handle:
        file_handle.next()
        return {'spot': float(next(file_handle))}


def parse_rm_duped_metrics_file_pe(rmDup_file):
    """
    Parses the rmdup file (tabbed file containing
     barcodes found ('randomer'),
     number of reads found ('total_count'),
     number of reads removed ('removed_count')

    :param rmDup_file: basestring
        filename of the rmDup file
    :return count_dict: dict
        dictionary containing sums of total, removed,
        and usable (total - removed)
    """
    # print("returning number of reads mapped in rmduped file")
    # return parse_rm_duped_metrics_file_se(rmDup_file)
    ########################################
    ### TODO: FIX THIS ###
    
    try:
        df = pd.read_csv(rmDup_file, sep="\t")
        return {
            "total_count": sum(df.total_count),
            "removed_count": sum(df.removed_count),
            "Usable reads": sum(df.total_count) - sum(df.removed_count)
        }
    except Exception as e:
        print(e)
        return {
            "total_count": None,
            "removed_count": None,
            "Usable reads": None
        }
    
    
def parse_rm_duped_metrics_file_se(rmDup_file, label="Usable reads"):
    """
    Parses the BAM file to return the number of reads mapped.
    (in the future when we umi tools produce stats pages this will change)

    :param rmDup_file: basestring
        BAM file name
    :return:
    """
    check_for_index(rmDup_file)
    samfile = pysam.AlignmentFile(rmDup_file)
    return {
        label: samfile.mapped
    }


def check_for_index(bamfile):

    """

    Checks to make sure a BAM file has an index, if the index does not exist it is created

    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file

    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))

    if os.path.exists(bamfile + ".bai"):
        return

    if not bamfile.endswith(".bam"):
        raise NameError("file %s not of correct type" % (bamfile))
    else:
        process = call(["samtools", "index", str(bamfile)])

        if process == -11:
            raise NameError("file %s not of correct type" % (bamfile))


def main():
    parser = argparse.ArgumentParser(description="Make a summary csv files of all eclash metrics")
    parser.add_argument("--analysis_dir", help="analysis directory", required=False, default="./")
    parser.add_argument("--output_csv", help="output csv filename", required=False, default="./eclashsummary.csv")
    parser.add_argument("--number_usable", help="number of usable reads", required=False, type=float, default=1000000)
    parser.add_argument("--percent_usable", help="percent of usable reads", required=False, type=float, default=0.7)
    parser.add_argument("--pipeline", help="either [total] or targeted", required=False, default='total')
    args = parser.parse_args()
    
    if args.pipeline == 'total':
        print("Aggregating results from total ChimClip run.")
        clashseq_metrics(
            args.analysis_dir,
            args.output_csv,
            args.percent_usable,
            args.number_usable,
            cutadapt_round1_suffix="*fqTr.metrics",
            repeat_mapped_suffix="*.filtered.chimeric_candidates.STARLog.final.out",
            repeat_mapped_eclip_suffix="*.umi.r1.fqTrTr.fq.sorted.STARLog.final.out",
            genome_mapped_suffix="*.sorted.repeat-unmapped.STARLog.final.out",
            genome_mapped_eclip_suffix="*.eclip.repeat-unmapped.STARLog.final.out",
            rm_dup_suffix="*.fq.sorted.genome-mappedSoSo.rmDupSo.bam",
            rm_dup_eclip_suffix="*.umi.r1.fqTrTr.fq.sorted.eclip.genome-mappedSoSo.rmDupSo.bam",
            peak_suffix=None,
        )
    elif args.pipeline == 'targeted':
        print("Aggregating results from targeted ChimClip run.")
        clashseq_metrics(
            args.analysis_dir,
            args.output_csv,
            args.percent_usable,
            args.number_usable,
            cutadapt_round1_suffix="*fqTr.metrics",
            repeat_mapped_suffix="*.primer.fq.sorted.STARLog.final.out",
            genome_mapped_suffix="*.repeat-unmapped.fq.sorted.STARLog.final.out",
            rm_dup_suffix="*.sorted.genome-mappedSoSo.rmDupSo.bam",
            peak_suffix="*.rmDupSo.peakClusters.merged.bed",
        )
if __name__ == '__main__':
    main()