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


def create_mock_transcript(name="TR1",
                           chrom="CHR1",
                           startPos=3,
                           cigar="8M7D6M2I2M11D7M",
                           direction="+"):
    """
    Utility function to create mock transcripts
    without the need to specify all inputs.
    """

    testTranscript = Transcript(name=name,
                                chrom=chrom,
                                startPos=startPos,
                                cigar=cigar,
                                direction=direction)

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

    testPos = create_mock_transcript(cigar="8M7D6M2I2M11D7M",
                                     startPos=3,
                                     direction="+")

    assert isinstance(testPos.conversionTree, IntervalTree)

    assert testPos.conversionTree == IntervalTree([Interval(0, 8, 3),
                                                   Interval(8, 14, 10),
                                                   Interval(14, 15, 9.1),
                                                   Interval(15, 16, 8.2),
                                                   Interval(16, 18, 8),
                                                   Interval(18, 25, 19)])

    testNeg = create_mock_transcript(cigar="8M7D6M2I2M11D7M",
                                     startPos=43,
                                     direction="-")

    assert testNeg.conversionTree == IntervalTree([Interval(0, 7, 43),
                                                   Interval(7, 9, 32),
                                                   Interval(9, 10, 33.1),
                                                   Interval(10, 11, 34.2),
                                                   Interval(11, 17, 34),
                                                   Interval(17, 25, 27)])


def test_translate_coordinates_pos_strand():
    """
    First Example
    COORD       0    5    10   15   20     25   30   35   40   45   50
    GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
    TR1         5' GTCATGTA-------CTAGCCGGTA-----------AGATAAT 3'
                   |    |           |    |               |   |
                   0    5           10   15              20  24
    CIGAR = 8M7D6M2I2M11D7M
    """

    testTranscript = create_mock_transcript(cigar="8M7D6M2I2M11D7M",
                                            direction="+")

    with pytest.raises(ValueError):
        testTranscript.translate_coordinates(-1)

    with pytest.raises(ValueError):
        testTranscript.translate_coordinates(25)

    expected_3 = {'name': 'TR1', 'inputPos': 3,
                  'chrom': 'CHR1', 'refPos': 6, 'direction': "+"}

    assert testTranscript.translate_coordinates(3) == expected_3

    expected_14 = {'name': 'TR1', 'inputPos': 14,
                   'chrom': 'CHR1', 'refPos': 23.1, 'direction': "+"}

    assert testTranscript.translate_coordinates(14) == expected_14

    expected_15 = {'name': 'TR1', 'inputPos': 15,
                   'chrom': 'CHR1', 'refPos': 23.2, 'direction': "+"}

    assert testTranscript.translate_coordinates(15) == expected_15

    expected_16 = {'name': 'TR1', 'inputPos': 16,
                   'chrom': 'CHR1', 'refPos': 24, 'direction': "+"}

    assert testTranscript.translate_coordinates(16) == expected_16

    """
    Second example
    COORD       0    5      10   15   20     25   30   35   40   45   50
    GENOME:CHR1 ACTGTCATGT--GTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
    TR1           5' CA---CAGTTT-------CTATAGCTAG-----------ATTTAGC 3'
                     |       |           |    |               |   |
                     0       5           10   15              20  24
    CIGAR = 2M3D2I4M7D1M2I7M11D7M
    """

    testTranscript2 = create_mock_transcript(startPos=5,
                                             cigar="2M3D2I4M7D1M2I7M11D7M",
                                             direction="+")

    expected2_3 = {'name': 'TR1', 'inputPos': 2,
                   'chrom': 'CHR1', 'refPos': 9.1, 'direction': "+"}

    assert testTranscript2.translate_coordinates(2) == expected2_3

    expected2_8 = {'name': 'TR1', 'inputPos': 8,
                   'chrom': 'CHR1', 'refPos': 21, 'direction': "+"}

    assert testTranscript2.translate_coordinates(8) == expected2_8

    expected2_24 = {'name': 'TR1', 'inputPos': 24,
                    'chrom': 'CHR1', 'refPos': 46, 'direction': "+"}

    assert testTranscript2.translate_coordinates(24) == expected2_24


def test_translate_coordinates_neg_strand():
    """
    First Example
    COORD       0    5    10   15   20     25   30   35   40   45   50
    GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATGTATTAGCTAG
    TR1         3' GTCATGTA-------CTAGCCGGTA-----------AGATGTA 5'
                   |   |           |    |               |    |
                   24  20          15   10              5    0
    CIGAR = 8M7D6M2I2M11D7M
    """
    testTranscript = create_mock_transcript(name="TR3",
                                            startPos=43,
                                            cigar="8M7D6M2I2M11D7M",
                                            direction="-")

    expected_0 = {'name': 'TR3', 'inputPos': 0,
                  'chrom': 'CHR1', 'refPos': 43, 'direction': "-"}

    assert testTranscript.translate_coordinates(0) == expected_0

    expected_9 = {'name': 'TR3', 'inputPos': 9,
                  'chrom': 'CHR1', 'refPos': 24.1, 'direction': "-"}

    assert testTranscript.translate_coordinates(9) == expected_9

    expected_10 = {'name': 'TR3', 'inputPos': 10,
                   'chrom': 'CHR1', 'refPos': 24.2, 'direction': "-"}

    assert testTranscript.translate_coordinates(10) == expected_10

    expected_24 = {'name': 'TR3', 'inputPos': 24,
                   'chrom': 'CHR1', 'refPos': 3, 'direction': "-"}

    assert testTranscript.translate_coordinates(24) == expected_24

    with pytest.raises(ValueError):
        testTranscript.translate_coordinates(-1)

    with pytest.raises(ValueError):
        testTranscript.translate_coordinates(25)

    """
    Second example
    COORD       0    5      10   15   20     25   30   35   40   45   50
    GENOME:CHR1 ACTGTCATGT--GTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
    TR1           3' CA---CAGTTT-------CTATAGCTAG-----------ATTTAGC 5'
                     |       |           |    |               |   |
                     24      19          14   9               4   0
    CIGAR = 2M3D2I4M7D1M2I7M11D7M
    """

    testTranscript2 = create_mock_transcript(startPos=46,
                                             cigar="2M3D2I4M7D1M2I7M11D7M",
                                             direction="-")

    expected2_0 = {'name': 'TR1', 'inputPos': 0,
                   'chrom': 'CHR1', 'refPos': 46, 'direction': "-"}

    assert testTranscript2.translate_coordinates(0) == expected2_0

    expected2_14 = {'name': 'TR1', 'inputPos': 14,
                    'chrom': 'CHR1', 'refPos': 22.1, 'direction': "-"}

    assert testTranscript2.translate_coordinates(14) == expected2_14

    expected2_22 = {'name': 'TR1', 'inputPos': 22,
                    'chrom': 'CHR1', 'refPos': 10.2, 'direction': "-"}

    assert testTranscript2.translate_coordinates(22) == expected2_22


def test_import_transcripts():

    testMapper = TranscriptMapper()

    expected = [{'name': 'TR1',
                 'chrom': 'CHR1',
                 'startPos': 3,
                 'cigar': '8M7D6M2I2M11D7M',
                 'direction': "+"},
                {'name': 'TR2',
                 'chrom': 'CHR2',
                 'startPos': 10,
                 'cigar': '20M',
                 'direction': "+"},
                {'name': 'TR3',
                 'chrom': 'CHR1',
                 'startPos': 43,
                 'cigar': '8M7D6M2I2M11D7M',
                 'direction': "-"}]

    result = testMapper.get_transcript_info_from_file(exampleTranscriptFile)

    for i in range(0, len(result)):
        assert result[i] == expected[i]

    with pytest.raises(FileNotFoundError):
        testMapper.get_transcript_info_from_file("./non_existing_file.tsv")


def test_import_query():

    testMapper = TranscriptMapper()

    expected = [{'name': 'TR1', 'queryPos': 4},
                {'name': 'TR2', 'queryPos': 0},
                {'name': 'TR3', 'queryPos': 0},
                {'name': 'TR1', 'queryPos': 13},
                {'name': 'TR2', 'queryPos': 10},
                {'name': 'TR3', 'queryPos': 9}]

    result = testMapper.get_query_from_file(exampleQueryFile)

    for i in range(0, len(result)):
        assert result[i] == expected[i]

    with pytest.raises(FileNotFoundError):
        testMapper.get_query_from_file("./non_existing_file.tsv")


def test_run_single_query():
    testMapper = TranscriptMapper()
    testMapper.import_transcripts(exampleTranscriptFile)

    singleQueryResult = testMapper.run_single_query("TR1", 24)

    expectedResult = {'name': 'TR1', 'inputPos': 24,
                      'chrom': 'CHR1', 'refPos': 43, 'direction': "+"}

    assert singleQueryResult == expectedResult


def test_check_transcript_line():

    validLine = "TR1\tCHR1\t3\t8M7D6M2I2M11D7M\t+\n"
    result1 = TranscriptMapper.check_transcript_line(validLine)
    expected1 = ["TR1", "CHR1", "3", "8M7D6M2I2M11D7M", "+"]

    assert result1 == expected1

    invalidLine1 = "TR1\tCHR1\t3\t8M7D6M2I2M11D7M\t+\t1\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_transcript_line(invalidLine1)

    invalidLine2 = "TR1\tCHR1\tB\t8M7D6M2I2M11D7M\t+\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_transcript_line(invalidLine2)

    invalidLine3 = "TR1\tCHR1\t1\t8M7D6M2I2M11D7M\t*\n"
    with pytest.raises(Exception):
        TranscriptMapper.check_transcript_line(invalidLine3)


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
