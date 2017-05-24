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

from parse_cutadapt import parse_cutadapt_file
from qcsummary_rnaseq import rnaseq_metrics_df, parse_star_file


def clipseq_metrics_csv(analysis_dir, output_csv):
    df = clipseq_metrics_df(analysis_dir=analysis_dir, iclip=True)
    df.to_csv(output_csv)


def clipseq_metrics_df(analysis_dir,
                       iclip=False,
                       num_seps=None,
                       sep=".",
                       percent_usable=.01,
                       number_usable=500000):
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
    if num_seps is None:
        num_seps = 2 if iclip else 1

    ###########################################################################
    # get file paths
    ################
    print("---")
    print("GET FILE PATHS FOR ECLIP")

    #cutadapt_round2_files = glob.glob(os.path.join(
    #    analysis_dir, "*.adapterTrim.round2.metrics"))
    cutadapt_round2_files = glob.glob(os.path.join(
        analysis_dir, "*fqTrTr.metrics"))
    print("cutadapt_round2_files:", cutadapt_round2_files)

    #rmRep_mapping_files = glob.glob(os.path.join(
    #    analysis_dir, "*.adapterTrim.round2.rep.bamLog.final.out"))
    rmRep_mapping_files = glob.glob(os.path.join(
        analysis_dir, "*fqTrTrU-SoMa.metrics"))
    print("rmRep_mapping_files:", rmRep_mapping_files)

    #rm_duped_files = glob.glob(os.path.join(analysis_dir, "*rmRep.rmDup.metrics"))
    rm_duped_files = glob.glob(os.path.join(analysis_dir, "*fqTrTrMa.metrics"))
    print("rm_duped_files:", rm_duped_files)
    # hack for new data
    #if len(rm_duped_files) == 0:
        #rm_duped_files = glob.glob(os.path.join(analysis_dir, "*r2.rmDup.metrics"))
    # hack for new data
    #if len(rm_duped_files) == 0:
        #rm_duped_files = glob.glob(os.path.join(analysis_dir, "*.rmDup.metrics"))

    #spot_files = glob.glob(os.path.join(
    #    analysis_dir, "*peaks.metrics"))
    spot_files = glob.glob(os.path.join(
        analysis_dir,"*fqTrTrU-SoMaCoSoV2Cl.bed.log"))
    print("spot_files:", spot_files)

    #peaks_files = glob.glob(os.path.join(
    #    analysis_dir, "*.peaks.bed"))
    peaks_files = glob.glob(os.path.join(
        analysis_dir, "*fqTrTrU-SoMaCoSoV2Cl.bed"))
    print("peaks_files:", peaks_files)

    print("---")
    ###########################################################################

    ###########################################################################
    # get file names
    ################
    cutadapt_round2_names = get_names(cutadapt_round2_files, num_seps, sep)
    rmRep_mapping_names = get_names(rmRep_mapping_files, num_seps, sep)
    rm_duped_names = get_names(rm_duped_files, num_seps, sep)
    spot_names = get_names(spot_files, num_seps, sep)
    peaks_names = get_names(peaks_files, num_seps, sep)
    ###########################################################################

    ###########################################################################
    # make dataframes
    #################
    cutadapt_round2_df = pd.DataFrame(
        {name: parse_cutadapt_file(cutadapt_file)
         for name, cutadapt_file in cutadapt_round2_names.items()}
    ).transpose()
    cutadapt_round2_df.columns = ["{} Round 2".format(col)
                                  for col in cutadapt_round2_df.columns]

    rmRep_mapping_df = pd.DataFrame(
        {name: parse_star_file(star_file)
         for name, star_file in rmRep_mapping_names.items()}
    ).transpose()
    rmRep_mapping_df.columns = ["{} rmRep".format(col)
                                for col in rmRep_mapping_df.columns]

    rm_duped_df = pd.DataFrame(
        {name: parse_rm_duped_metrics_file(rm_duped_file)
         for name, rm_duped_file in rm_duped_names.items()}
    ).transpose()

    # FIXME TEMPORARILY COMMENTED OUT
    # spot_df = pd.DataFrame(
    #     {name: parse_peak_metrics(spot_file)
    #      for name, spot_file in spot_names.items()}
    # ).transpose()

    peaks_df = pd.DataFrame(
        {name: {"Num Peaks": len(pybedtools.BedTool(peaks_file))}
         for name, peaks_file in peaks_names.items()}
    ).transpose()
    ###########################################################################

    ###########################################################################
    # get rnaseq metrics dataframe
    ##############################
    combined_df = rnaseq_metrics_df(analysis_dir, num_seps, sep)
    ###########################################################################

    ###########################################################################
    # merge dataframes
    ##################
    combined_df = pd.merge(combined_df, cutadapt_round2_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rm_duped_df,
                           left_index=True, right_index=True, how="outer")
    # FIXME TEMPORARILY COMMENTED OUT
    # combined_df = pd.merge(combined_df, spot_df,
    #                        left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, peaks_df,
                           left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rmRep_mapping_df,
                           left_index=True, right_index=True, how="outer")
    ###########################################################################


    ###########################################################################
    # compute useful stats
    ######################
    combined_df['Uniquely Mapped Reads'] = \
        combined_df['Uniquely Mapped Reads'].astype(float)
    try:
        combined_df["Percent Usable / Input"] = \
            (combined_df['Usable Reads'] / combined_df['Uniquely Mapped Reads'])
        combined_df["Percent Usable / Mapped"] = \
            (combined_df['Usable Reads'] / combined_df['Input Reads'])
        combined_df['Passed QC'] = \
            (combined_df['Usable Reads'] > number_usable) & \
            (combined_df['Percent Usable / Mapped'] > percent_usable)
    except ZeroDivisionError:
        print("passing on ZeroDivisionError")
        pass
    ###########################################################################

    return combined_df


def get_names(files, num_seps, sep):
    ################################
    """
    Given a list of files return that files base name and the path to that file
    Args:
        files: list of files
        num_seps: int number of seperators in to call the real name
        sep: str seperator to split on
    Returns: dict basename to file
    """
    dict_basename_to_file = {sep.join(os.path.basename(file).split(sep)[0: num_seps]): file
                             for file in files}
    return dict_basename_to_file


def parse_peak_metrics(fn):
    #######################
    with open(fn) as file_handle:
        file_handle.next()
        return {'spot': float(file_handle.next())}


def parse_rm_duped_metrics_file(rmDup_file):
    ########################################
    try:
        df = pd.read_csv(rmDup_file, sep="\t")
        return {
            "total_count": sum(df.total_count),
            "removed_count": sum(df.removed_count),
            "Usable Reads": sum(df.total_count) - sum(df.removed_count)}
    except Exception as e:
        print(e)
        return {
            "total_count": None,
            "removed_count": None,
            "Usable Reads": None
        }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Make a summary csv files of all eclip metrics")
    parser.add_argument("--analysis_dir", help="analysis directory", required=False, default="./")
    parser.add_argument("--output_csv", help="output csf filename", required=False, default="./eclipqcsummary.csv")
    args = parser.parse_args()
    # print("args:", args)
    clipseq_metrics_csv(args.analysis_dir, args.output_csv)
