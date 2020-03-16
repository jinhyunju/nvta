# Transcript Mapping Utilities 

`nvta` is a python library for mapping transcript coordinates to the reference genome based on its CIGAR string.

# Installation

Using an environment that has `python3`, execute the following command to install the package.

`pip install -e /<path>/<to>/<package>`

If you extracted the contents of the tar ball you can point to the directory, or you can also use the tar ball directly.

`pip install nvta-0.1.tar.gz`

# External Dependencies

This package uses the following dependencies:

- `intervaltree`
- `pytest`

# Tutorial

For a detailed tutorial of the features, please use this jupyter notebook.

`./tutorial/transcript_utils_tutorial.ipynb`

# Quick start example

## Processing multiple transcripts and queries from files

1. Import `TranscriptMapper` and set example inputs

```python
from nvta.transcript_utils import TranscriptMapper

fileTranscriptInput = "./tests/resources/example_transcript_input.tsv"
fileQueryInput = "./tests/resources/example_query.tsv"
outputFile = "./test_result.tsv"
```

2. Initiate a `TranscriptMapper` object and import information from files.

```python
transcriptMapper = TranscriptMapper()

# import transcripts and queries
transcriptMapper.import_transcripts(fileTranscriptInput)
transcriptMapper.import_queries(fileQueryInput)
```

3. Execute all queries and export results to output file.

```python
# execute all queries
transcriptMapper.run_all_queries(showResults=True)

# write results to file
transcriptMapper.export_query_results(outputFile=outputFile)
```

## Example usage of the `Transcript` class with direct user input

```python
from nvta.transcript_utils import Transcript

# construct a Transcript object
singleTranscript = Transcript(name="TR1",
                              startPos=3,
                              chrom="CHR1",
                              cigar="8M7D6M2I2M11D7M",
                              direction="+")

# map a transcript coordinate to the reference
singleTranscript.translate_coordinates(4)
```

# Original Problem Statement

### 5' to 3' mapping direction (`+`)
The objective is to write software that translates transcript coordinates to genomic coordinates. For example, consider the simple transcript TR1, which aligns to the genome as follows:

```
COORD       0    5    10   15   20     25   30   35   40   45   50
GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
TR1         5' GTCATGTA-------CTAGCCGGTA-----------AGATAAT 3'
               |    |           |    |               |   |
               0    5           10   15              20  24
```

We can compactly express this alignment in the same way that we compactly represent a read alignment in the SAM/BAM format : using a position and CIGAR string. In this case, the (0-based) position is CHR1:3, and the CIGAR string is `8M7D6M2I2M11D7M`. For this exercise, you may assume that the transcript is always mapped from genomic 5’ to 3’.

The objective is then to translate a (0-based) transcript coordinate to a (0 based) genome coordinate. For example the fifth base in TR1 (i.e. TR1:4) maps to genome coordinate CHR1:7. Similarly, TR1:13 maps to CHR1:23 and TR1:14 maps to an insertion immediately after CHR1:23, which is represented as CHR1:23.1 in this implementation.

### 3' to 5' mapping direction (`-`)

For the same transcript that maps in the opposite direction the CIGAR string will remain the same as `8M7D6M2I2M11D7M`, but the transcript coordinates are going to be reversed as shown below.

```
COORD       0    5    10   15   20     25   30   35   40   45   50
GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATGTATTAGCTAG
TR1         3' GTCATGTA-------CTAGCCGGTA-----------AGATGTA 5'
               |   |           |    |               |    |
               24  20          15   10              5    0
```

Here we would map the start of the transcript to CHR1:43. The ninth base would map to the first base of the insertion after CHR1:24, and will be represented as CHR1:24.1 whereas the tenth base will map to the second base of the insertion as CHR1:24.2. The end of the transcript, 24, will map to CHR1:3.


# Implementation Details

## Inputs and Outputs

### Inputs
`Transcript input`: A five column (tab-separated) file containing the transcripts. 
- The first column is the transcript name (`string`)
- The second column shows the chromosome that the transcript maps to (`string`)
- The third column shows the 0-based starting position of the transcript (5' end) on the reference. (`int`)
- The fourth column has the CIGAR string indicating the mapping of the transcript to the reference. (`string`)
- The fifth column should have `+` for a 5'-3' direction and a `-` for a 3'-5' direction transcript. (`string`)

Example :

```
TR1	CHR1	3	8M7D6M2I2M11D7M	+
```

---

`Query input`: A two column (tab-separated) file indicating a set of queries. 
- The first column shows the name of the transcript. (`string`)
- The second column shows the 0-based transcript coordinate that should be translated to a reference position. (`int`)

Example:

```
TR1	4
```

---

### Outputs

`Query results`: The output is a five column tab separated file with one row for each of the query results. 
- The first two columns are exactly the two columns from the query input file 
- The third column shows the chromosome name (`string`)
- The fourth column shows the reference base the query mapped to (`int` for regular, `float` for insertions)
- The fifth column shows the direction of the transcript (`string`)

Example:

```
TR1	4	CHR1	7	+
```

---

### Additional detail for handling insertions 

To easily distinguish insertion bases, which do not have a direct mapping to the reference, I have implemented a modified output for the insertion coordinates. The output coordinate will have the format of `<ref_base_before_insertion>.<number_of_inserted_base>`. For example, the first position of the insertion (14) in the problem statement would map to `23.1` since the reference position before the start of the insertion is 23 and it is the 1st base in the insertion. The 1-based coordinate for insertion was specifically chosen to avoid confusions such as `23` vs`23.0` if a 0-based insertion coordinate was chosen. 


## Key Assumptions

- Transcript CIGAR string does not contain S (soft-clipped), H (hard-clipped) operations. This is assuming that the provided transcripts are not directly derived from read mappings, but have been processed to start and end at bases actually mapping to the reference.

- Transcript CIGAR string does not contain P (padding) operations. This is assuming that multiple sequence alignments are not required for this task.

- The user is aware of the size of the reference and is responsible to check whether the query results are going to generate out of bounds values on the reference given a transcript and cigar. 

- CIGAR strings are always represented in the 5'-3' direction of the reference. 

## Strengths and Weaknesses

- This solution supports the use cases described in the problem statement of reading in transcripts and queries from a file, and also provides the user the flexibility to create custom workflows by having an independent `Transcript` class. This also makes it easier to create a parallelized solution for dealing with a large number of unique transcripts.

- Information is stored in objects that are mainly default data types with the exception of `IntervalTree`s making it easier for the user to retrieve and inspect the content.

- Solution is fast for cases where there are long stretches of the same operation in a CIGAR string by grouping them into a single interval that translates the coordinates by a simple addition operation. 

- Since all of the intervals for a transcript are mutually exclusive, alternatives to `IntervalTree` might yield better performance. However, the use of `IntervalTree` makes it easier to inspect the translation since it has both the intervals and the values in one place. 

- Performance might be sub-optimal for cases where there are a lot of small stretches of different CIGAR operations.

- When the goal is to support a lot of queries that are not necessarily unique, pre-computing the translation and saving it in a look-up table might yield better performance.


## Testing

- Unit tests covering the functionality of the implementation can be found in `./tests/test_transcript_utils.py`.

- Both the functionality and error reporting are covered in the unit tests.

- All tests can can be executed with the following command: 

```
pytest tests
```

## Additional Thoughts

- To produce real-world data using external data sources, one might download cDNA sequences from public sources such as Ensembl (using their FTP download features) and map them to a reference in house to generate the CIGAR strings. This would allow the user to not be tied to a specific version of the reference. 
- For dealing with very long CIGAR strings, one could split the transcript between long stretches of deletions (likely introns) to reduce the amount of data that is associated to a single transcript. 
- To handle a large number of CIGAR strings, the data could be partitioned into chromosomes and sub regions within the chromosome, so that each subset is at a reasonable size. 
- If one expects a large number of queries that are not unique, pre-computing the translation of coordinates and having a look up table might be more efficient (as mentioned in the strengths/weaknesses section). One could also save the result of a specific query to a look up table once it has been executed and directly look it up from the table when it is queried again. Having values pre-computed for the most frequently queried genes/transcripts only and doing the translation on the fly for other less frequently queried genes/transcripts could be a happy medium as well.
