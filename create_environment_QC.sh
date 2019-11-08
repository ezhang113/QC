#!/bin/bash

conda remove -y --name qc --all;

conda create -y --name qc python=3.6;
source activate qc;

conda install -c bioconda -c anaconda \
matplotlib pandas numpy seaborn bedtools pybedtools pysam xlsxwriter
