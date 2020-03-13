import re
import numpy as np
import logging
from intervaltree import IntervalTree

logger = logging.getLogger(__name__)


class Transcript():
    """
    Transcript class with the following functionalities
    
    - Process a CIGAR string from a single transcript
    - Generate an Intervaltree to transform coordinates
    - Transform coordinates to the reference based on input
    
    """
    
    def __init__(self, name, chrom, startPos, cigar):
        """
        Initiate transcript object with following inputs
        
        :param name: string describing the transcript name
        :param chrom: string for transcript chromosome (ex. chr1, chr2)
        :param startPos: int specifying start position of transcript on reference
        :param cigar: string containing the CIGAR string of the transcript 
        :returns: none
        :rtype: none
        """
        self.name = name
        self.chrom = chrom
        self.startPos = startPos
        # validate cigar
        self.cigar = self.verify_cigar(cigar)
        # create interval tree based on cigar
        self.conversionTree = self.process_cigar()
        # get max transcript position
        self.transcriptEnd = self.conversionTree.end() - 1 
        
    def get_info(self):
        """
        Accessor method to get the transcript information

        :returns: A dictionary with the transcript 
                  name, chromosome, start position, and CIGAR string
        :rtype: dict
        """
        transcriptInfo = {'name' : self.name, 
                           'chrom' : self.chrom, 
                           'startPos' : self.startPos, 
                           'cigar' : self.cigar}
        return transcriptInfo

    @staticmethod
    def verify_cigar(cigar):
        """
        Method to verify a given CIGAR string

        :param cigar: string containing the CIGAR string of the transcript 
        :returns: Input cigar string if it passed all checks
        :rtype: string
        """

        validCigarChar = ["M","I","D","N","S","H","P","X","="]
        matches = re.findall(r'([0-9]+)([A-Z|a-z|:/?#@!$&'"'"'()*+,;=%[]{1})', cigar)
        logger.debug("Processing CIGAR string {}".format(cigar))

        reconstruct = ""
        for cigarEntry in matches:
            opInt = int(cigarEntry[0])
            opChar = cigarEntry[1]

            if opChar not in validCigarChar:
                # constant adjustment
                logger.error("Invalid CIGAR string character detected {}".format(opChar))
                raise ValueError("Invalid CIGAR string character")

            # need to check if H and S
            if opChar in ['H', 'S', 'P']:
                logger.error("CIGAR parser does not support logic for {}".format(opChar))
                raise ValueError("Unsupported CIGAR character")
            reconstruct += cigarEntry[0] + cigarEntry[1]

        if len(reconstruct) != len(cigar):
            logger.error("Parsed information does not match original CIGAR, potentially malformatted CIGAR string.")
            logger.error("Input CIGAR = {}".format(cigar))
            logger.error("Parsed CIGAR = {}".format(reconstruct))
            raise ValueError("Malformatted CIGAR string")
        
        return cigar
    
    def process_cigar(self):
        """
        Process a CIGAR string and create IntervalTree for coordinate translation
        
        :returns: IntervalTree containing intervals of the transcript coordinates 
                  that can be translated to the reference coordinates by adding 
                  the value saved for the interval
        :rtype: IntervalTree
        """
        logger.info("Start Processing of Transcript CIGAR")

        refTransform = self.startPos
        cigar = self.cigar
        
        tStart = 0
        tEnd = 0
        conversionTree = IntervalTree()

        matches = re.findall(r'([0-9]+)([MIDNSHPX=]{1})', cigar)
        
        for cigarEntry in matches:
            logger.debug("Processing CIGAR entry {}".format(cigarEntry))

            opInt = int(cigarEntry[0])
            opChar = cigarEntry[1]

            if opChar in ['M','=','X']:
                # constant adjustment
                tEnd = tStart + opInt
                conversionTree[tStart:tEnd] = refTransform

            elif opChar in ['N','D']:
                # increase adjustment by gap size 
                refTransform = refTransform + opInt

            elif opChar in ['I']:
                # handle insertion to reference
                for i in range(0,opInt):
                    
                    tEnd = tStart + 1
                    refTransform = refTransform - 1
                    conversionTree[tStart:tEnd] = float(str(refTransform)+"."+str(i + 1))
                    tStart = tEnd

            elif opChar in ['H','S','P']:
                continue

            else:
                # catch unknown
                logger.error("Encountered invalid character in CIGAR {}".format(opChar))

            tStart = tEnd
        logger.info("Finished Processing of Transcript CIGAR")
        return conversionTree


    def translate_coordinates(self, inputPosition):
        """
        Translate a transcript position to the reference position
        
        :param inputPosition: int specifying the transcript coordinate to be translated
        :returns: dictionary containing the transcript name (name)
                  input position (inputPos), chromosome (chrom), and translated ref position (refPos)
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
            adjustedCoordinate = inputPosition + posAdjustment
            
        results = {'name' : self.name, 
                   'inputPos' : inputPosition, 
                   'chrom' : self.chrom, 
                   'refPos' : adjustedCoordinate}
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
        self.queryInfo = []

    def get_queries(self):
        """
        Accessor for query info

        :returns: list of dictionaries containing imported queries
        :rtype: list
        """
        return self.queryInfo

    def get_transcripts(self):
        """
        Accessor for transcripts

        :returns: dictionary of imported transcripts
        :rtype: dict
        """
        return self.queryInfo

    def _build_transcripts(self, inputFile):
        tmpTxDict = {}
        transcriptInfo = self.get_transcript_info_from_file(inputFile)
        for item in transcriptInfo :
            tmpTxDict[item['name']] = Transcript(**item)
        return tmpTxDict
    
    def import_transcripts(self, inputFile):
        self.transcripts = self._build_transcripts(inputFile)
    
    def import_queries(self, inputFile):
        self.queryInfo = self.get_query_from_file(inputFile)

    def run_all_queries(self, showResults = False):
        results = [self.run_single_query(**x) for x in self.queryInfo]
        self.queryResults = results
        if showResults:
            return self.queryResults

    def run_single_query(self, name, queryPos):
        try:
            result = self.transcripts[name].translate_coordinates(queryPos)
        except KeyError:
            logger.error("Input query contains a transcript that has not been loaded")
            raise ValueError
        return result
    
    def export_query_results(self, outputFile):
        # check if file exists
        with open(outputFile, 'w') as f:
            for item in self.queryResults:
                f.write("{name}\t{inputPos}\t{chrom}\t{refPos}\n".format(**item))

    @staticmethod
    def get_transcript_info_from_file(inputFile):
        inputTranscripts = []
        logger.info("Reading Transcript Information from {}".format(inputFile))
        with open(inputFile, 'r') as f:
            for line in f:
                # check if 4 tab separated entries are present
                # and whether the 3rd entry can be converted to int
                lineSplit = TranscriptMapper._check_transcript_line(line)
                
                transcriptInfo = {'name' : lineSplit[0], 
                                  'chrom' : lineSplit[1],
                                  'startPos' : int(lineSplit[2]),
                                  'cigar' : lineSplit[3]}
                inputTranscripts.append(transcriptInfo)
        return inputTranscripts
    
    @staticmethod
    def get_query_from_file(inputFile):
        inputQuery = []
        logger.info("Reading Query Information from {}".format(inputFile))
        with open(inputFile, 'r') as f:
            for line in f:
                lineSplit = TranscriptMapper._check_query_line(line)
                inputQuery.append({'name' : lineSplit[0], 
                                   'queryPos' : int(lineSplit[1])})
        return inputQuery
    
    @staticmethod
    def _check_transcript_line(line):
        lineSplit = line.strip().split("\t")
        
        if len(lineSplit) != 4:
            raise Exception("Input file must be tab separated and have 4 entries.")
        
        try:
            int(lineSplit[2])
        except:
            raise Exception("Third column needs to be an integer")

        return lineSplit

    @staticmethod
    def _check_query_line(line):
        lineSplit = line.strip().split("\t")
        
        if len(lineSplit) != 2:
            raise Exception("Input file must be tab separated and have 2 entries.")
        
        try:
            int(lineSplit[1])
        except:
            raise Exception("Second entry needs to be an integer")
        return lineSplit
