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
# from qcsummary_rnaseq import rnaseq_metrics_df, parse_star_file
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
Percent usable / mapped
Percent Usable / Input""".split("\n")

clash_slim_qc_metrics = [
    "Initial reads num", "Reads after cutadapt 1", "STAR genome input reads",
    "STAR genome uniquely mapped", "STAR genome uniquely mapped %",
    'Number of reads mapped to too many loci',
    '% of reads unmapped: too short', '% of reads mapped to too many loci', "Usable reads",
    "Percent usable / mapped", "Percent Usable / Input"
]

def write_clipseq_metrics_csv(df, output_csv):
    """
    Writes CSV file with clipseq metrics.
    Writes "annotated" (abridged with the most useful QC stats) file.
    """
    df.to_csv(output_csv)
    df[clash_slim_qc_metrics].to_csv(os.path.splitext(output_csv)[0] + ".annotated.csv")
    
    
def clashseq_metrics(analysis_dir, output_csv, percent_usable, number_usable):
    
    # TODO: remove iclip param when new nomenclature is finalized.
    df = clipseq_metrics_df(
        analysis_dir=analysis_dir,
        percent_usable=percent_usable,
        number_usable=number_usable,
        iclip=False,
    )
    df = df[CLASH_SE_ORDER]
    
    write_clipseq_metrics_csv(df, output_csv)
    

def clipseq_metrics_df(
        analysis_dir, percent_usable,
        number_usable,
        iclip=False, num_seps=None,
        sep=".",
        cutadapt_round1_suffix="*fqTr.metrics",
        cutadapt_round2_suffix="*fqTrTr.metrics",
        genome_mapped_suffix="*.filtered.chimeric_candidates.STARLog.final.out",
        rm_dup_suffix="*.filtered.chimeric_candidates.STARAligned.outSo.rmDup.bam",
        # peak_suffix="*.peakClusters.bed",
    ):
    #######################################
    """
    Reports all clip-seq metrics in a given analysis directory
    outputs must follow gabes naming clipseq pipeline / naming conventions"
    Args:
        analysis_dir:
        iclip:
        num_seps:
        sep:
        percent_usable:
        number_usable:
    Returns:
    """
    # TODO: fix prefix name separator
    if num_seps is None:
        num_seps = 3 if iclip else 3

    cutadapt_round1_names, \
    cutadapt_round2_names, \
    genome_mapped_names, \
    rm_duped_names = get_all_names(
        analysis_dir=analysis_dir,
        cutadapt_round1_suffix=cutadapt_round1_suffix,
        cutadapt_round2_suffix=cutadapt_round2_suffix,
        genome_mapped_suffix=genome_mapped_suffix, 
        rm_dup_suffix=rm_dup_suffix,
        # peak_suffix=peak_suffix,
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
    if len(cutadapt_round2_names) > 0:
        cutadapt_round2_df = pd.DataFrame(
            {
                name: pc.parse_cutadapt_file(cutadapt_file, paired_end=False)
                for name, cutadapt_file in cutadapt_round2_names.items()
            }
        ).transpose()

        cutadapt_round2_df.columns = [
            "{} Round 2".format(col) for col in cutadapt_round2_df.columns
        ]
    if len(genome_mapped_suffix) > 0:
        star_df = pd.DataFrame(
            {
                name: parse_star_file(genome_mapped_file) 
                for name, genome_mapped_file in genome_mapped_names.items()
            }
        ).transpose()
    if len(rm_duped_names) > 0:
        rm_duped_df = pd.DataFrame(
            {
                name: parse_rm_duped_metrics_file_se(rm_duped_file)
                for name, rm_duped_file in rm_duped_names.items()
            }
        ).transpose()
    else:
        rm_duped_df = pd.DataFrame(index=['Usable reads', 'removed_count', 'total_count']).T
    # if len(peaks_names) > 0:
    #     peaks_df = pd.DataFrame(
    #         {
    #             name: {"Clipper peaks num": len(pybedtools.BedTool(peaks_file))}
    #             for name, peaks_file in peaks_names.items()
    #         }
    #     ).transpose()
        
    ###########################################################################

    ###########################################################################
    # get the base dataframe rnaseq metrics dataframe
    ##############################
    ###########################################################################

    ###########################################################################
    # merge dataframes
    ##################
    # combined_df = pd.merge(combined_df, cutadapt_round2_df,
    #                        left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, star_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rm_duped_df,
                           left_index=True, right_index=True, how="outer")
    # combined_df = pd.merge(combined_df, peaks_df,
    #                        left_index=True, right_index=True, how="outer")
    
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
            "Number of input reads": "STAR genome input reads",
            "STAR genome uniquely mapped number": "STAR genome uniquely mapped",
            "Reads Written Round 1": "Reads after cutadapt 1"
        }
    )
    # compute useful stats
    ######################
    combined_df['STAR genome uniquely mapped'] = combined_df['STAR genome uniquely mapped'].astype(float)
    combined_df['Initial reads num'] = combined_df['Initial reads num'].astype(float)
    try:
        combined_df["Percent usable / mapped"] = (
                    combined_df['Usable reads'] / combined_df['STAR genome uniquely mapped'])
        combined_df["Percent Usable / Input"] = (combined_df['Usable reads'] / combined_df['Initial reads num'])
        combined_df["Percent Repetitive"] = 1 - (
                    combined_df['STAR genome input reads'] / combined_df['Reads after cutadapt 1'].astype(float))
        combined_df["Repetitive Reads"] = combined_df['Reads after cutadapt 1'] - combined_df['STAR genome input reads']
        # combined_df['Passed basic QC'] = (combined_df['Usable reads'] > number_usable) & (
        #             combined_df['Percent usable / mapped'] > percent_usable)

    except ZeroDivisionError:
        print("passing on ZeroDivisionError")
        pass

    return combined_df


def get_all_names(
        analysis_dir,
        cutadapt_round1_suffix,
        cutadapt_round2_suffix,
        genome_mapped_suffix, 
        rm_dup_suffix,
        # peak_suffix,
        sep,
        num_seps
    ):
    ###########################################################################
    # get file paths
    ################

    cutadapt_round1_files = glob.glob(os.path.join(analysis_dir, cutadapt_round1_suffix))
    cutadapt_round2_files = glob.glob(os.path.join(analysis_dir, cutadapt_round2_suffix))
    genome_mapped_files = glob.glob(os.path.join(analysis_dir, genome_mapped_suffix))
    rm_duped_files = glob.glob(os.path.join(analysis_dir, rm_dup_suffix))
    # peaks_files = glob.glob(os.path.join(analysis_dir, peak_suffix))
    ###########################################################################

    ###########################################################################
    # get file names
    ################
    cutadapt_round1_names = get_names(cutadapt_round1_files, num_seps, sep)
    cutadapt_round2_names = get_names(cutadapt_round2_files, num_seps, sep)
    genome_mapped_names = get_names(genome_mapped_files, num_seps, sep)
    rm_duped_names = get_names(rm_duped_files, num_seps, sep)
    # peaks_names = get_names(peaks_files, num_seps, sep)
    ###########################################################################
    return cutadapt_round1_names, cutadapt_round2_names, genome_mapped_names, rm_duped_names # , peaks_names


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
    
    
def parse_rm_duped_metrics_file_se(rmDup_file):
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
        "Usable reads": samfile.mapped
    }


def build_second_mapped_from_master(df):
    second_mapped = df[[
        '% of reads unmapped: too short',
        '% of reads mapped to too many loci',
        '% of reads unmapped: too many mismatches',
        'STAR genome uniquely mapped %',
        'Percent usable / mapped'
    ]].fillna('0')
    for col in second_mapped.columns:
        try:
            second_mapped[col] = second_mapped[col].apply(
                lambda x: float(x.strip('%')) / 100
            )
        except AttributeError:
            second_mapped[col] = second_mapped[col].astype(float)
    return second_mapped


def build_peak_df_from_master(df):
    peaks = df[[
        'Clipper peaks num',
    ]]

    return peaks


def build_raw_number_from_master(df):
    num = df[[
        'Usable reads',
        'STAR genome input reads',
        'STAR genome uniquely mapped',
        'Repetitive Reads'
    ]]
    return num


def parse_star_file(star_file_name):
    with open(star_file_name) as star_file:
        star_dict = {}
        star_dict["Started job on"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Started mapping on"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Finished on"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Mapping speed, Million of reads per hour"] = next(star_file).strip().split("|")[1].strip()
        next(star_file)
        star_dict["Number of input reads"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Average input read length"] = float(next(star_file).strip().split("|")[1].strip())
        next(star_file)
        star_dict["STAR genome uniquely mapped number"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["STAR genome uniquely mapped %"] = float(next(star_file).strip().split("|")[1].strip()[:-1])
        star_dict["Average mapped length"] = float(next(star_file).strip().split("|")[1].strip())
        star_dict["Number of splices: Total"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Number of splices: Annotated (sjdb)"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Number of splices: GT/AG"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Number of splices: GC/AG"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Number of splices: AT/AC"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Number of splices: Non-canonical"] = int(next(star_file).strip().split("|")[1].strip())
        star_dict["Mismatch rate per base, percent"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Deletion rate per base"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Deletion average length"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Insertion rate per base"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Insertion average length"] = next(star_file).strip().split("|")[1].strip()
        next(star_file)
        star_dict["Number of reads mapped to multiple loci"] = next(star_file).strip().split("|")[1].strip()
        star_dict["% of reads mapped to multiple loci"] = next(star_file).strip().split("|")[1].strip()
        star_dict["Number of reads mapped to too many loci"] = next(star_file).strip().split("|")[1].strip()
        star_dict["% of reads mapped to too many loci"] = next(star_file).strip().split("|")[1].strip()
        next(star_file)
        star_dict["% of reads unmapped: too many mismatches"] = next(star_file).strip().split("|")[1].strip()
        star_dict["% of reads unmapped: too short"] = next(star_file).strip().split("|")[1].strip()
        star_dict["% of reads unmapped: other"] = next(star_file).strip().split("|")[1].strip()
        
    return star_dict

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
    args = parser.parse_args()
    
    clashseq_metrics(
        args.analysis_dir,
        args.output_csv,
        args.percent_usable,
        args.number_usable,
    )

if __name__ == '__main__':
    main()