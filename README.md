# Transcript Mapping Utilities 

`nvta` is a python library for mapping transcript coordinates to the reference genome based on its CIGAR string.

# To Do

[v] Need to add strand to input
[v] Need to add test for reverse
[v] Need to add strand to output
[] Need to update readme with direction

# Installation

`pip install -e /<path>/<to>/<package>`

# External Dependencies

This package uses the following dependencies:

- `intervaltree`
- `pytest`

# Quick start example

## Processing multiple transcripts and queries from files

```python
from nvta.transcript_utils import TranscriptMapper

fileTranscriptInput = "./tests/resources/example_transcript_input.tsv"
fileQueryInput = "./tests/resources/example_query.tsv"
outputFile = "./test_result.tsv"

transcriptMapper = TranscriptMapper()

# import transcripts and queries
transcriptMapper.import_transcripts(fileTranscriptInput)
transcriptMapper.import_queries(fileQueryInput)

# execute all queries
transcriptMapper.run_all_queries(showResults = True)

# write results to file
transcriptMapper.export_query_results(outputFile= outputFile)
```

## Use single transcript with direct user input

```python
from nvta.transcript_utils import Transcript

singleTranscript = Transcript(name="TR1",
                              startPos=3,
                              chrom="CHR1",
                              cigar="8M7D6M2I2M11D7M")

singleTranscript.translate_coordinates(4)

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

## Inputs and Outputs

### Inputs
- Transcript input: A four column (tab-separated) file containing the transcripts. The first column is the transcript name, and the remaining three columns indicate it’s genomic mapping: chromosome name, 0-based starting position on the chromosome, and CIGAR string indicating the mapping.

- Query input: A two column (tab-separated) file indicating a set of queries. The first column is a transcript
name, and the second column is a 0-based transcript coordinate.

### Outputs

- Query results: The output is a four column tab separated file with one row for each of the input queries. The first two columns are exactly the two columns from the second input file, and the remaining two columns are the chromosome name and chromosome coordinate, respectively.

- Additional detail for handling insertions: To easily distinguish insertions that do not have a direct mapping to the reference from other bases, I have implemented a modified output for the insertion coordinates. The output coordinate will have the format of `<ref_base_before_insertion>.<number_of_inserted_base>`. For example, the first position of the insertion (14) in the problem statement would map to `23.1` since the reference position before the start of the insertion is 23 and it is the 1st base in the insertion. The 1-based coordinate for insertion was specifically chosen to avoid confusions such as `23` vs`23.0` if a 0-based insertion coordinate was chosen. 


## Key Assumptions

- Transcript CIGAR string does not contain S (soft-clipped), H (hard-clipped) operations. This is assuming that the transcripts provided are not directly derived from read mappings, but have been processed to start and end at bases actually mapping to the reference. 

- Transcript CIGAR string does not contain P (padding) operations. This is assuming that multiple sequence alignments are not required for this task.

- The user is aware of the size of the reference and is responsible to check whether the supplied transcript and cigar is not going to generate out of bounds values for a given coordinate translation. 

## Strengths and Weaknesses

- This solution supports the use cases described in the problem statement of reading in transcripts and queries from a file, and also provides the user the flexibility to create custom workflows by having an independent `Transcript` class. It also makes it easier to create a parallelized solution for dealing with a large number of unique transcripts.

- Information is stored in objects that are mainly default data types with the exception of `IntervalTree`s making it easier for the user to retrieve and inspect the content.

- Solution is fast for cases where there are long stretches of the same operation in a CIGAR string by grouping them into a single interval that translates the coordinates by a simple addition operation. 

- Since all of the intervals for a transcript are mutually exclusive, alternatives to `IntervalTree` might yield better performance. However, the use of `IntervalTree` makes it easier to inspect the translation since it has both the intervals and the values in one place. 

- Performance might be sub-optimal for cases where there are a lot of small stretches of different CIGAR operations.

- When the goal is to support a lot of queries that are not neccessary unique, pre-computing the translation and saving it in a look-up table might yield better performance.


## Testing

- Unit tests covering the functionality of the implementation can be found under `/tests/` and can be executed with the following command `pytest tests`.


## Reverse Mapping


```
COORD       0    5    10   15   20     25   30   35   40   45   50
GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
TR1            GTCATGTA-------CTAGCCGGTA-----------AGATAAT
               |   |           |    |               |    |
               24  20          15   10              5    0
```

0 - 43  (+ 43)
6 - 37  (+ 31)
11D
43 - 11 - 7 = 32 - 7 = 25
7 - 25  (+ 18)

9 - 26.1
10 - 26.2

We can compactly express this alignment in the same way that we compactly represent a read alignment in the SAM/BAM format : using a position and CIGAR string. In this case, the (0-based) position is CHR1:3, and the CIGAR string is 8M7D6M2I2M11D7M . For this exercise, you may assume that the transcript is always mapped from genomic 5’ to 3’.

The objective is then to translate a (0-based) transcript coordinate to a (0 based) genome coordinate. For example the fifth base in TR1 (i.e. TR1:4) maps to genome coordinate CHR1:7. Similarly, TR1:13 maps to CHR1:23 and TR1:14 maps to an insertion immediately before CHR1:24.
