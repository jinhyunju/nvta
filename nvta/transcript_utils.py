import re
import os
import logging
from intervaltree import IntervalTree

logger = logging.getLogger(__name__)


class Transcript():
    """
    Transcript class with the following functionalities

    - Process a CIGAR string from a single transcript
    - Generate an Intervaltree to transform coordinates
    - Transform coordinates to the reference based on input
    - Supporting both 5'-3' and 3'-5' mapping.

    """

    def __init__(self, name, chrom, startPos, cigar, direction="+"):
        """
        Initiate transcript object with following inputs

        :param name: string describing the transcript name
        :param chrom: string for transcript chromosome (ex. chr1, chr2)
        :param startPos: int specifying start position on reference
        :param cigar: string containing the CIGAR string of the transcript
        :param direction: string specifying the transcript direction.
                          '+' for 5'-3' and '-' for 3'-5' directions.
                          default is '+'
        :return: none
        :rtype: none
        """
        self.name = name
        self.chrom = chrom
        self.startPos = startPos
        self.direction = direction
        # validate cigar
        self.cigar = self.verify_cigar(cigar)
        # create interval tree based on cigar
        self.conversionTree = self.process_cigar()
        # get max transcript position
        self.transcriptEnd = self.conversionTree.end() - 1

    def get_info(self):
        """
        Accessor method to get the transcript information

        :return: A dictionary with the transcript
                 name, chromosome, start position, CIGAR string, and direction
        :rtype: dict
        """
        transcriptInfo = {
                          'name': self.name,
                          'chrom': self.chrom,
                          'startPos': self.startPos,
                          'cigar': self.cigar,
                          'direction': self.direction
                          }
        return transcriptInfo

    @staticmethod
    def verify_cigar(cigar):
        """
        Method to verify a given CIGAR string

        :param cigar: string containing the CIGAR string of the transcript
        :return: Input cigar string if it passed all checks
        :rtype: string
        """

        validCigarChar = ["M", "I", "D", "N", "S", "H", "P", "X", "="]
        matches = re.findall(r'([0-9]+)([A-Z|a-z|:/?#@!$&'"'"'()*+,;=%[]{1})',
                             cigar)
        logger.debug("Processing CIGAR string {}".format(cigar))

        reconstruct = ""
        for cigarEntry in matches:
            opInt = int(cigarEntry[0])
            opChar = cigarEntry[1]

            if opChar not in validCigarChar:
                # constant adjustment
                logger.error("""Invalid CIGAR string character
                                detected {}""".format(opChar))
                raise ValueError("Invalid CIGAR string character")

            # need to check if H and S
            if opChar in ['H', 'S', 'P']:
                logger.error("""CIGAR parser does not support
                                logic for {}""".format(opChar))
                raise ValueError("Unsupported CIGAR character")
            reconstruct += cigarEntry[0] + cigarEntry[1]

        if len(reconstruct) != len(cigar):
            logger.error("""Parsed information does not match original CIGAR,
                            potentially malformatted CIGAR string.""")
            logger.error("Input CIGAR = {}".format(cigar))
            logger.error("Parsed CIGAR = {}".format(reconstruct))
            raise ValueError("Malformatted CIGAR string")

        return cigar

    def process_cigar(self):
        """
        Process a CIGAR string and create an IntervalTree
        for coordinate translation

        :return: IntervalTree containing intervals of
                 the transcript coordinates that can be
                 translated to the reference coordinates
                 when adjusted by the value saved for the interval
        :rtype: IntervalTree
        """
        logger.info("Start Processing of Transcript CIGAR")

        refTransform = self.startPos
        cigar = self.cigar
        direction = self.direction

        tStart = 0
        tEnd = 0
        conversionTree = IntervalTree()
        matches = re.findall(r'([0-9]+)([MIDNX=]{1})', cigar)

        if direction == "-":
            # if direction is 3'-5' reverse the cigar operation order
            matches = reversed(matches)

        for cigarEntry in matches:
            logger.debug("Processing CIGAR entry {}".format(cigarEntry))

            opInt = int(cigarEntry[0])
            opChar = cigarEntry[1]

            if opChar in ['M', '=', 'X']:
                # constant adjustment for matching bases
                tEnd = tStart + opInt
                conversionTree[tStart:tEnd] = refTransform

            elif opChar in ['N', 'D']:
                # increase adjustment to account for gaps
                if direction == "-":
                    # if 3'-5' reduce correction factor
                    refTransform = refTransform - opInt
                else:
                    # if 5'-3' increase correction factor
                    refTransform = refTransform + opInt

            elif opChar in ['I']:
                # special logic to handle insertion to reference
                for i in range(0, opInt):

                    tEnd = tStart + 1
                    if direction == "-":
                        #  if 3'-5' add insertion base back to correction
                        refTransform = refTransform + 1
                    else:
                        # if 5'-3' subtract insertion base from correction
                        refTransform = refTransform - 1
                    conversionTree[tStart:tEnd] = float(str(refTransform) +
                                                        "." +
                                                        str(i + 1))
                    tStart = tEnd

            else:
                # catch unknown if it wasn't caught in CIGAR validation
                logger.error("""Encountered invalid character
                                in CIGAR {}""".format(opChar))
                raise ValueError("Unknown CIGAR operation {}".format(opChar))
            # advance to next transcript interval
            tStart = tEnd

        logger.info("Finished Processing of Transcript CIGAR")

        return conversionTree

    def translate_coordinates(self, inputPosition):
        """
        Translate a transcript position to the reference position

        :param inputPosition: int specifying the transcript
                              coordinate to be translated
        :return: dictionary containing the transcript name (name)
                 input position (inputPos), chromosome (chrom),
                 translated ref position (refPos),
                 and direction (direction).
        :rtype: dict
        """
        if inputPosition < 0:
            logger.error("Please use a valid position")
            raise ValueError("Negative position given.")

        if inputPosition > self.transcriptEnd:
            logger.error("Input position out of transcript bounds.")
            raise ValueError("Position exceeding transript length.")

        intervalSet = self.conversionTree[inputPosition]
        if len(intervalSet) > 1:
            logger.error("Overlapping intervals detected!")
            raise Exception("Error in interval tree")

        elif len(intervalSet) == 0:
            logger.error("No interval for position detected!")
            raise Exception("Error in interval tree")

        else:
            posAdjustment = intervalSet.pop().data

        if self.direction == "-":
            adjustedCoordinate = posAdjustment - inputPosition
        else:
            adjustedCoordinate = inputPosition + posAdjustment

        # this is not super clean, but was required to deal with
        # python's floating point arithmetic limitations
        if float(adjustedCoordinate).is_integer():
            refCoordinate = int(adjustedCoordinate)
        else:
            refCoordinate = float(format(adjustedCoordinate, '.1f'))

        results = {'name': self.name,
                   'inputPos': inputPosition,
                   'chrom': self.chrom,
                   'refPos': refCoordinate,
                   'direction': self.direction}

        return results


class TranscriptMapper():
    """
    Transcript mapper class that functions as a utility wrapper for Transcripts
    Features include:

    - Import single/multiple transcript information from a file
    - Import single/multiple queries from a file
    - Execute all queries
    - Write query results to a file
    """

    def __init__(self):
        self.transcripts = {}
        self.queryResults = []
        self.queries = []

    def get_queries(self, index=None):
        """
        Accessor for query info

        :param index: int specifying the 0-based index of a single query
        :return: list of dictionaries containing imported queries
        :rtype: list
        """
        if index is None:
            # if no index is supplied return all queries
            return self.queries
        else:
            if index >= 0 & index < len(self.queries):
                return self.queries[index]
            else:
                raise Exception("Query index out of bounds")

    def get_transcripts(self, name=None):
        """
        Accessor for transcripts

        :param name: string specifying a transcript name to retrieve
        :return: dictionary of imported transcripts
        :rtype: dict
        """
        if name is None:
            # if no name is supplied return all transcripts
            return self.transcripts
        else:
            if name in self.transcripts:
                return self.transcripts[name]
            else:
                raise Exception("Transcript not found.")

    def _build_transcripts(self, inputFile):
        """
        Internal wrapper for retrieving transcript information.

        :param inputFile: string containing path to inputFile.
        :return: dict of Transcripts with transcript names as keys
        :rtype: dict
        """
        tmpTxDict = {}
        transcriptInfo = self.get_transcript_info_from_file(inputFile)

        for item in transcriptInfo:
            tmpTxDict[item['name']] = Transcript(**item)
        return tmpTxDict

    def import_transcripts(self, inputFile):
        """
        Main method for importing transcripts from files.

        :param inputFile: string containing path to input file.
        :return: none
        """
        self.transcripts = self._build_transcripts(inputFile)

    def import_queries(self, inputFile):
        """
        Main method for importing queris from files.

        :param inputFile: string containing path to input file.
        :return: none
        """
        self.queries = self.get_query_from_file(inputFile)

    def run_all_queries(self, showResults=False):
        """
        Method to run all queries imported to object.

        :param showResults: bool if set to True return query results
        :return: list of dictionaries containing query results
        :rtype: list
        """
        results = [self.run_single_query(**x) for x in self.queries]
        self.queryResults = results
        if showResults:
            return results

    def run_single_query(self, name, queryPos):
        """
        Method to run a single query

        :param name: string specifying the transcript name
        :param queryPos: int specifying the transcript position to translate.
        :return: single query result
        :rtype: dict
        """
        try:
            result = self.transcripts[name].translate_coordinates(queryPos)
        except KeyError:
            logger.error("""Input query contains a transcript
                            that has not been loaded""")
            raise ValueError
        return result

    def export_query_results(self, outputFile):
        """
        Method to exort saved query results to file.

        :param outputFile: string specifying the output file
                           (path to output file must exist)
        :return: none
        """
        # check if file exists
        with open(outputFile, 'w') as f:
            for item in self.queryResults:
                f.write("""{name}\t{inputPos}\t{chrom}\t
                           {refPos}\t{direction}\n""".format(**item))

    @staticmethod
    def get_transcript_info_from_file(inputFile):
        """
        Static method to read transcript information from a file

        :param inputFile: string containing path to input file.
        :return: list of dictionaries containing transcript information
        :rtype: list
        """
        inputTranscripts = []
        if not os.path.exists(inputFile):
            raise FileNotFoundError("{} not found".format(inputFile))

        logger.info("Reading Transcript Information from {}".format(inputFile))
        with open(inputFile, 'r') as f:
            for line in f:
                # check if 4 tab separated entries are present
                # and whether the 3rd entry can be converted to int
                lineSplit = TranscriptMapper.check_transcript_line(line)

                transcriptInfo = {'name': lineSplit[0],
                                  'chrom': lineSplit[1],
                                  'startPos': int(lineSplit[2]),
                                  'cigar': lineSplit[3],
                                  'direction': lineSplit[4]}
                inputTranscripts.append(transcriptInfo)
        return inputTranscripts

    @staticmethod
    def get_query_from_file(inputFile):
        """
        Static method to read query information from a file

        :param inputFile: string containing path to input file.
        :return: list of dictionaries containing query information
        :rtype: list
        """
        if not os.path.exists(inputFile):
            raise FileNotFoundError("{} not found".format(inputFile))

        inputQuery = []
        logger.info("Reading Query Information from {}".format(inputFile))
        with open(inputFile, 'r') as f:
            for line in f:
                lineSplit = TranscriptMapper.check_query_line(line)
                inputQuery.append({'name': lineSplit[0],
                                   'queryPos': int(lineSplit[1])})
        return inputQuery

    @staticmethod
    def check_transcript_line(line):
        """
        Helper method to verify a line of a transcript input file

        :param line: string containing a single line from a transcript file.
        :return: list of elements in line separated by tab.
        :rtype: list
        """
        lineSplit = line.strip().split("\t")

        if len(lineSplit) != 5:
            raise Exception("Input file must have 5 tab separated entries.")

        try:
            int(lineSplit[2])
        except ValueError:
            raise Exception("Third column needs to be an integer")

        if lineSplit[4] not in ["+", "-"]:
            raise Exception("""Fifth column needs to be '+' or '-'
                               indicating the direction of the transcript""")

        return lineSplit

    @staticmethod
    def check_query_line(line):
        """
        Helper method to verify a line of a query input file

        :param line: string containing a single line from a query file.
        :return: list of elements in line separated by tab.
        :rtype: list
        """
        lineSplit = line.strip().split("\t")

        if len(lineSplit) != 2:
            raise Exception("Input file must have 2 tab separated entries.")

        try:
            int(lineSplit[1])
        except ValueError:
            raise Exception("Second entry needs to be an integer")
        return lineSplit
