# `nvta` package with `transcript_utils`

A python library for mapping transcript coordinates to the reference genome based on its CIGAR string.

# Installation

`pip install -e /<path>/<to>/<package>`

# Quick start example

```
from nvta.transcript_utils import Transcript, TranscriptMapper

fileTranscriptInput = "./tests/resources/example_transcript_input.tsv"
fileQueryInput = "./tests/resources/example_query.tsv"
outputFile = "./test_result.tsv"

transcriptMapper = TranscriptMapper()

transcriptMapper.import_transcripts(fileTranscriptInput)

transcriptMapper.import_queries(fileQueryInput)

transcriptMapper.run_all_queries(showResults = True)

transcriptMapper.export_query_results(outputFile= outputFile)
```

# Implementation Details

## Original Problem Statment

The objective is to write software that translates transcript coordinates to genomic coordinates. For example, consider the simple transcript TR1, which aligns to the genome as follows:

```
COORD       0    5    10   15   20     25   30   35   40   45   50
GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
TR1            GTCATGTA-------CTAGCCGGTA-----------AGATAAT
               |    |           |    |               |   |
               0    5           10   15              20  24
```
We can compactly express this alignment in the same way that we compactly represent a read alignment in the SAM/BAM format : using a position and CIGAR string. In this case, the (0-based) position is CHR1:3, and the CIGAR string is 8M7D6M2I2M11D7M . For this exercise, you may assume that the transcript is always mapped from genomic 5’ to 3’.

The objective is then to translate a (0-based) transcript coordinate to a (0 based) genome coordinate. For example the fifth base in TR1 (i.e. TR1:4) maps to genome coordinate CHR1:7. Similarly, TR1:13 maps to CHR1:23 and TR1:14 maps to an insertion immediately before CHR1:24.


## Key Assumptions

- Transcript CIGAR string does not contain S (soft-clipped), H (hard-clipped) operations. This is assuming that the transcripts provided are not directly derived from read mappings, but have been processed to start at a base which actually maps to the reference. 

- Transcript CIGAR string does not contain P (padding) operations. This is assuming that multiple sequence alignments are not required for this task.

- 
## Strengths and Weaknesses
- Solution is intuitive.
- Speed up can be achieved since transcript coordinates intervals are not overlapping, so more efficient data structures could be used in place of the IntervalTree.