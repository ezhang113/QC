#!/usr/bin/env bash

cd results/intermediates
qcsummaryeclip
mv eclipqcsummary.csv ../qc-summary.csv
cd -
