#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool


baseCommand: [qcsummary_eclip.py]


arguments: [
  --bw_pos,
  $(inputs.bam.nameroot).posbw,
  --bw_neg,
  $(inputs.bam.nameroot).negbw
  ]

inputs:

  bam:
     type: File
     format: http://edamontology.org/format_2572
     inputBinding:
       position: 1
       prefix: --bam
     #secondaryFiles:
     #  - ".bai"

  bai:
    type: File
    format: http://edamontology.org/format_3327
    inputBinding:
      position: 2
      prefix: --bai

  chromsizes:
    type: File
    inputBinding:
      position: 3
      prefix: --genome

outputs:

  posbw:
    type: File
    format: http://edamontology.org/format_3006
    outputBinding:
      glob: $(inputs.bam.nameroot).posbw

  negbw:
    type: File
    format: http://edamontology.org/format_3006
    outputBinding:
      glob: $(inputs.bam.nameroot).negbw
