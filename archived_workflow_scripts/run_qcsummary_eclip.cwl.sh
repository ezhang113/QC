#!/usr/bin/env cwl-runner

export PATH=${PWD}:$PATH
echo $PATH
cwl-runner qcsummary_eclip.cwl qcsummary_job.yaml
