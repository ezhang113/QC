# QC summary parser

This script is designed to run after the completion of an eCLIP pipeline run.
It will search an analysis folder for metrics and summary files, and
combine all into a single comma-separated file for easier parsing.

Usage:

```
qcsummary_eclip.py \
--analysis_dir results/   # folder containing all eCLIP pipeline intermediates \
--output_csv qcsummary.csv  # output filename \
--number_usable 1500000   # number of usable reads (typically ~1.5M)\
--percent_usable 0.5   # percent of PCR-deduped / total mapped reads \
--peak_threshold 1000  #
```

### VERY GENERAL GUIDELINES
For human (hg19) datasets, first-pass QC requires 1.5M reads and 1000 peaks.
This is dependent on the RBP binding behavior 