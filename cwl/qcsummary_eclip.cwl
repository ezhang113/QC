#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool


baseCommand: [qcsummary_eclip.py]


arguments: [
  --output_csv,
  $(inputs.analysis_dir.)/eclipsummary.csv
  ]

inputs:

  analysis_dir:
     type: Directory
     inputBinding:
       position: 1
       prefix: --analysis_dir

  number_usable:
    type: int
    inputBinding:
      position: 2
      prefix: --number_usable

  percent_usable:
    type: float
    inputBinding:
      position: 3
      prefix: --percent_usable

  peak_threshold:
    type: int
    inputBinding:
      position: 4
      prefix: --peak_threshold

outputs:

  summary_csv:
    type: File
    outputBinding:
      glob: $(inputs.analysis_dir)/eclipsummary.csv

  summary_figure:
    type: File
    outputBinding:
      glob: $(inputs.analysis_dir)/eclipsummary.png
