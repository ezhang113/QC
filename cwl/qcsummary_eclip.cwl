#!/usr/bin/env cwltool

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

baseCommand: [qcsummary_eclip.py]


arguments: [
  --output_csv,
  $(inputs.analysis_dir.path)/eclipsummary.csv
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

  output_csv:
    default: ""
    type: string
    inputBinding:
      position: 5
      prefix: --output_csv
      valueFrom: |
        ${
          if (inputs.output_csv == "") {
            return "qcsummary.csv";
          }
          else {
            return inputs.output_csv;
          }
        }


outputs:

  summary_csv:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_csv == "") {
            return "qcsummary.csv";
          }
          else {
            return inputs.output_csv;
          }
        }

  summary_figure:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_csv == "") {
            return "qcsummary.png";
          }
          else {
            return inputs.output_csv;
          }
        }
