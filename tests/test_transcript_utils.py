import os
import pytest
from intervaltree import IntervalTree, Interval
from nvta.transcript_utils import Transcript, TranscriptMapper

"""
Resources for testing
"""
resourceDir = "./tests/resources"

exampleTranscriptFile = os.path.join(resourceDir, 
                                     "example_transcript_input.tsv")

exampleQueryFile = os.path.join(resourceDir, 
                                "example_query.tsv")


"""
Template for creating a mock transcript
"""
def create_mock_transcript(name = "TR1", 
                           chrom = "CHR1", 
                           startPos = 3, 
                           cigar = "8M7D6M2I2M11D7M"):

    testTranscript = Transcript(name = name, 
                                chrom = chrom, 
                                startPos = startPos,
                                cigar = cigar)
    return testTranscript

"""
Testing starts here
"""

def test_verify_cigar():

    testCigar1 = "8M7D6M2I2M11D7M"
    resultCigar1 = Transcript.verify_cigar(testCigar1)

    assert resultCigar1 == testCigar1

    # test case for invalid character in cigar
    invalidCharCigar = "8M7D6M2I2M11D7M9O"
    with pytest.raises(ValueError):
        Transcript.verify_cigar(invalidCharCigar)

    malformatedCigar = "8MD7I"
    with pytest.raises(ValueError):
        Transcript.verify_cigar(malformatedCigar)


def test_process_cigar():

    testTranscript = create_mock_transcript(cigar = "8M7D6M2I2M11D7M")
    
    assert isinstance(testTranscript.conversionTree, IntervalTree)

    assert testTranscript.conversionTree[0] == {Interval(0, 8, 3)}


def test_translate_coordinates():
    
    testTranscript = create_mock_transcript(cigar = "8M7D6M2I2M11D7M")

    expected_3 = {'name': 'TR1', 'inputPos': 3, 'chrom': 'CHR1', 'refPos': 6}
    assert testTranscript.translate_coordinates(3) == expected_3

    expected_14 = {'name': 'TR1', 'inputPos': 14, 'chrom': 'CHR1', 'refPos': 23.1}
    assert testTranscript.translate_coordinates(14) == expected_14
    
    expected_15 = {'name': 'TR1', 'inputPos': 15, 'chrom': 'CHR1', 'refPos': 23.2}
    assert testTranscript.translate_coordinates(15) == expected_15

    expected_16 = {'name': 'TR1', 'inputPos': 16, 'chrom': 'CHR1', 'refPos': 24}
    assert testTranscript.translate_coordinates(16) == expected_16

    with pytest.raises(ValueError):
        testTranscript.translate_coordinates(-1)

    with pytest.raises(ValueError):
        testTranscript.translate_coordinates(25)

def test_import_transcripts():

    testMapper = TranscriptMapper()

    expected = [{'name' : 'TR1', 
                 'chrom' : 'CHR1',
                 'startPos' : 3,
                 'cigar' : '8M7D6M2I2M11D7M'},
                 {'name' : 'TR2', 
                 'chrom' : 'CHR2',
                 'startPos' : 10,
                 'cigar' : '20M'}]

    result = testMapper.get_transcript_info_from_file(exampleTranscriptFile)
    
    for i in range(0,len(result)):
        assert result[i] == expected[i]


def test_import_query():

    testMapper = TranscriptMapper()

    expected = [{'name' : 'TR1', 
                 'queryPos' : 4},
                 {'name' : 'TR2', 
                 'queryPos' : 0},
                 {'name' : 'TR1', 
                 'queryPos' : 13},
                 {'name' : 'TR2', 
                 'queryPos' : 10}
                 ]

    result = testMapper.get_query_from_file(exampleQueryFile)
    
    for i in range(0,len(result)):
        assert result[i] == expected[i]

def test_run_single_query():
    testMapper = TranscriptMapper()
    testMapper.import_transcripts(exampleTranscriptFile)
    
    singleQueryResult = testMapper.run_single_query("TR1", 24)

    expectedResult = {'name': 'TR1', 'inputPos': 24, 'chrom': 'CHR1', 'refPos': 42}

def test_check_transcript_line():

    validLine = "TR1\tCHR1\t3\t8M7D6M2I2M11D7M\n"
    result1 = TranscriptMapper.check_transcript_line(validLine)
    expected1 = ["TR1", "CHR1", "3", "8M7D6M2I2M11D7M"]

    assert result1 == expected1

    invalidLine1 = "TR1\tCHR1\t3\t8M7D6M2I2M11D7M\t1\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_transcript_line(invalidLine1)

    invalidLine2 = "TR1\tCHR1\tB\t8M7D6M2I2M11D7M\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_transcript_line(invalidLine2)

def test_check_query_line():

    validLine = "TR1\t10\n"
    result1 = TranscriptMapper.check_query_line(validLine)
    expected1 = ["TR1", "10"]

    assert result1 == expected1

    invalidLine1 = "TR1\t10\t8X\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_query_line(invalidLine1)

    invalidLine2 = "TR1\tU\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_query_line(invalidLine2)
