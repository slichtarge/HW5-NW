# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    sub_matrix_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10
    gap_extend = -1
    nw = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    score, seq1_align, seq2_align = nw.align(seq1, seq2)

    #shape
    assert nw._align_matrix.shape == (len(seq1) + 1, len(seq2) + 1)
    assert nw._gapA_matrix.shape == (len(seq1) + 1, len(seq2) + 1)
    assert nw._gapB_matrix.shape == (len(seq1) + 1, len(seq2) + 1)

    #base case
    assert nw._align_matrix[0, 0] == 0

    #first gap penalties
    assert nw._gapA_matrix[1, 0] == -11
    assert nw._gapB_matrix[0, 1] == -11

    #that the score is the actual max of all 3 matrices' last cell
    last_row, last_col = len(seq1), len(seq2)
    expected_max = max(nw._align_matrix[last_row, last_col], 
                        nw._gapA_matrix[last_row, last_col], 
                        nw._gapB_matrix[last_row, last_col])
    assert score == expected_max
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    sub_matrix_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10
    gap_extend = -1
    nw = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)
    score, seq3_align, seq4_align = nw.align(seq3, seq4)

        #alignments same length
    assert len(seq3_align) == len(seq4_align)

    #removing gaps should give us the oriignal sequences
    assert seq3_align.replace("-", "") == seq3
    assert seq4_align.replace("-", "") == seq4

    #we get expected result?
    assert seq3_align == "MAVHQLIRRP"
    assert seq4_align == "M---QLIRHP"





