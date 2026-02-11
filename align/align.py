# Importing Dependencies
from typing import Tuple
import numpy as np

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        


        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        
        rows = len(self._seqA) + 1
        cols = len(self._seqB) + 1

        # setting alignment and gap matrices to -infinity
        self._align_matrix = np.full((rows, cols), -np.inf) # M (alignment)
        self._gapA_matrix = np.full((rows, cols), -np.inf)  # gapA
        self._gapB_matrix = np.full((rows, cols), -np.inf)  # gapB

        # backtrace matrices (to store which matrix/direction we came from)
        # I'm defining -1 as blank, 0 as diagonal (alignment), 1 as up (gapA), and 2 as left (gapB)
        self._back = np.full((rows, cols), -1, dtype=int)
        self._back_A = np.full((rows, cols), -1, dtype=int)
        self._back_B = np.full((rows, cols), -1,  dtype=int)

        # base cases:
        self._align_matrix[0][0] = 0 #starting point

        # vertical gaps at start (seqA against gaps)
        for i in range(1, rows):        
            self._gapA_matrix[i][0] = self.gap_open + (i * self.gap_extend) #first gap is open + extend
            self._align_matrix[i][0] = self._gapA_matrix[i][0] #rest are extend
            self._back[i][0] = 1 

        # horizontal gaps at start (seqB against gaps)
        for j in range(1, cols):     
            self._gapB_matrix[0][j] = self.gap_open + (j * self.gap_extend)
            self._align_matrix[0][j] = self._gapB_matrix[0][j]
            self._back[0][j] = 2

        
        # TODO: Implement global alignment here

        #for every cell in ixj
        for i in range(1, rows):
            for j in range(1, cols):
                #gapA
                gapA_options = [
                    self._align_matrix[i-1][j] + self.gap_open + self.gap_extend, 
                    self._gapA_matrix[i-1][j] + self.gap_extend
                ]
                gapA_score = max(gapA_options)                
                self._gapA_matrix[i][j] = gapA_score
                self._back_A[i][j] = np.argmax(gapA_options) #keep track of which we chose (0 for M, 1 for gapA)

                #gapB
                gapB_options = [
                    self._align_matrix[i][j-1] + self.gap_open + self.gap_extend,
                    self._gapB_matrix[i][j-1] + self.gap_extend
                ]
                gapB_score = max(gapB_options)
                self._gapB_matrix[i][j] = gapB_score
                self._back_B[i][j] = np.argmax(gapB_options)*2  #added the *2 so that we get 0 for M, 2 for gapB

                #alignment
                similarity = self.sub_dict[(self._seqA[i-1], self._seqB[j-1])]

                M_options = [
                    self._align_matrix[i-1][j-1] + similarity,
                    self._gapA_matrix[i-1][j-1] + similarity,
                    self._gapB_matrix[i-1][j-1] + similarity
                ]
                M_score = max(M_options)
                self._align_matrix[i][j] = M_score
                self._back[i][j] = np.argmax(M_options)     #0 for M, 1 for gapA, 2 for gapB
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        i = len(self._seqA) #rows
        j = len(self._seqB) #cols

        #determine which matrix has the max score at the end
        final_score_options = [
            self._align_matrix[i][j],
            self._gapA_matrix[i][j],
            self._gapB_matrix[i][j]
        ]
        current_matrix = np.argmax(final_score_options)
        self.alignment_score = final_score_options[current_matrix] #the final alignment score is bottom right value of the max matrix

        seqA_align = []
        seqB_align = []

        while i > 0 or j > 0:

            #edge handling:
            if i == 0:          # if we hit the top row, we must move LEFT (aka gap in seqA)
                seqA_align.append('-')
                seqB_align.append(self._seqB[j-1])
                j -= 1
            elif j == 0:        # if we hit the left column, we must move UP (aka gap in seqB)
                seqA_align.append(self._seqA[i-1])
                seqB_align.append('-')
                i -= 1

            #normal backtrace
            else:
                #if align matrix, diagonal move
                if current_matrix == 0:
                    seqA_align.append(self._seqA[i-1])
                    seqB_align.append(self._seqB[j-1])
                    current_matrix = self._back[i][j]
                    i-=1
                    j-=1
                
                #if gapA matrix, move up (gap in seqB)
                elif current_matrix == 1:
                    seqA_align.append(self._seqA[i-1])
                    seqB_align.append('-')
                    current_matrix = self._back_A[i][j]
                    i-=1
                    

                #if gapB matrix, move to the left (gap in seqA)
                elif current_matrix == 2:
                    seqA_align.append('-')
                    seqB_align.append(self._seqB[j-1])
                    current_matrix = self._back_B[i][j]
                    j-=1

        #reverse and convert into string!
        self.seqA_align = "".join(reversed(seqA_align))
        self.seqB_align = "".join(reversed(seqB_align))

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header