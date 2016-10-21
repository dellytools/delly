#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "Delly"
label: "Delly tool"
cwlVersion: cwl:draft-3
description: |
  Tool:	     Delly
  Sources:   https://github.com/tobiasrausch/delly
  Summary:   Structural variant discovery by integrated paired-end and split-read analysis
requirements:
  - class: DockerRequirement
    dockerPull: "trausch/delly"
inputs:
  - id: "#type"
    type: string
    description: 'SV type'
    inputBinding:
      position: 1
      prefix: "-t"
  - id: "#genome"
    type: File
    description: 'Input genome file'
    inputBinding:
      position: 2
      prefix: "-g"
  - id: "#exclude"
    type: File
    description: 'Exclude file'
    inputBinding:
      position: 3
      prefix: "-x"
  - id: "#output"
    type: string
    description: 'Output file'
    inputBinding:
      position: 4
      prefix: "-o"
  - id: "#bams"
    type:
      type: array
      items: File
      secondaryFiles:
        - .bai
    inputBinding:
      position: 5
outputs: []
baseCommand: ["delly", "call"]
