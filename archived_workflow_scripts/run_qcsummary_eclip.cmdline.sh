#!/bin/bash

python qcsummary_eclip.py \
--peak_threshold 100 \
--number_usable 1000 \
--percent_usable 0.5 \
--analysis_dir /home/bay001/projects/codebase/QC/data/eclip_results-0.1.5/kbp550/results/intermediates \
--output_csv /home/bay001/projects/codebase/QC/data/eclip_results-0.1.5/kbp550/results/intermediates/newQC.csv
