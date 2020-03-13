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

