# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # # TODO Align all species to humans and print species in order of most similar to human BRD
    # # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    sub_matrix_file = "substitution_matrices/BLOSUM62.mat"
    gap_open = -10
    gap_extend = -1
    nw = NeedlemanWunsch(sub_matrix_file, gap_open, gap_extend)

    non_hs_seqs = [gg_seq, mm_seq, br_seq, tt_seq]
    non_hs_seq_labels = ["Gallus_gallus", "Mus_musculus", "Balaeniceps_rex", "tursiops_truncatus"]
    scores = []
    hs_seq_aligns = []
    non_hs_seq_aligns = []


    for non_hs_seq, non_hs_header in zip(non_hs_seqs, non_hs_seq_labels):
        score, hs_seq_align, non_hs_seq_align = nw.align(hs_seq, non_hs_seq)
        scores.append(score)
        hs_seq_aligns.append(hs_seq_align)
        non_hs_seq_aligns.append(non_hs_seq_align)

    sorted_lists = sorted(zip(scores, hs_seq_aligns, non_hs_seq_aligns, non_hs_seq_labels), key=lambda x: x[0], reverse=True) #zip and sort
    scores, hs_seq_aligns, non_hs_seq_aligns, non_hs_seq_labels = zip(*sorted_lists) #unzip
    scores, hs_seq_aligns, non_hs_seq_aligns, non_hs_seq_labels = list(scores), list(hs_seq_aligns), list(non_hs_seq_aligns), list(non_hs_seq_labels) #turn back to lists

    print("\nSpecies in order of most similar to human BRD: \n")
    for non_hs_seq_label in non_hs_seq_labels:
        print(non_hs_seq_label)
    print()

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    for score in scores:
        print(f"SCORE: {score}")
    print()
    

if __name__ == "__main__":
    main()
