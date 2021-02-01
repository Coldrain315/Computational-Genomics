import sys
from compsci260lib import *
from aligner_helpers import *
def run_ag_aligner_plus():
    """
    Align atpa_Hs.fasta and atpaEc.fasta and report the optimal alignment 
    score, the top-most alignment, and the bottom-most alignment.
    """

    #
    # YOUR CODE GOES HERE
    Ec = get_fasta_dict('atpa_Ec.fasta')[
        'sp|P0ABB0|ATPA_ECOLI ATP synthase subunit alpha OS=Escherichia coli (strain K12) OX=83333 GN=atpA PE=1 SV=1']
    Hs = get_fasta_dict('atpa_Hs.fasta')[
        'sp|P25705|ATPA_HUMAN ATP synthase subunit alpha, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1A PE=1 SV=1']
    #
    # seq1 = "GACTA"
    # seq2 = "GTTGAGTAC"
    match = 2
    mismatch = -1
    gap_penalty = 1
    affine_penalty = 11
    seq_type = validate_sequences(Ec, Hs)
    if seq_type == 1:
        # Both the sequences are DNA sequences so use the scores for match and
        # mismatch
        subst_matrix = create_subst_matrix_dna(match, mismatch)
        # print(subst_matrix)
    elif seq_type == 2:
        # Both the sequences are protein sequences so read in the BLOSUM62
        # substitution matrix
        subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    else:
        sys.exit('Input sequences are of different types: not both DNA or both protein')

    # Obtain a dictionary of scores for aligning a pair of characters
    subst_dict = create_subst_matrix_dict(subst_matrix)
    output_tuple = solve_ag_aligner_plus(Hs, Ec, subst_dict, gap_penalty, affine_penalty)
    print("The optimal alignment score is %d." % output_tuple[0])
    print("The top-most alignment sequence is %s and %s" % (output_tuple[1][0], output_tuple[1][1]))
    print("The bottom-most alignment sequence is %s and %s" % (output_tuple[2][0], output_tuple[2][1]))


def solve_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty):
    """The procedure for collecting the inputs, running the aligner,
    and returning the score and optimal alignments

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary representation
            of the substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character); this
            value should be positive because we will subtract it
        affine_penalty (int): affine penalty; as a positive integer

    Returns a tuple of:
        (the optimal alignment score as an int,
         the top-most alignment achieving this score as a tuple of strings,
         the bottom-most alignment achieving this score as a tuple of strings)

        Example output:

            (6, ("AT-AGG", "ATCCGG"), ("ATA-GG", "ATCCGG"))
    """

    #
    # YOUR CODE GOES HERE
    list1 = []
    # initialize two two-D array with n rows. and n is the length of seq1 + 1
    mid_table = [list(list1) for i in range(len(seq1) + 1)]
    upper_table = [list(list1) for i in range(len(seq1) + 1)]
    left_table = [list(list1) for i in range(len(seq1) + 1)]
    pointer_table = [list(list1) for i in range(len(seq1) + 1)]
    # go through all positions of two sequences. Plus one because there is an extra list and a column for '-'
    for i in range(len(seq1) + 1):
        mid_row = []  # initialize empty rows for dp table
        upper_row = []
        left_row = []
        pointer_row = [list(list1) for i in range(len(seq2) + 1)]  # initialize an array for pointer table
        for j in range(len(seq2) + 1):
            if i > 0 and j > 0:
                # for the F table
                left = max(left_table[i][j - 1] - gap_penalty, mid_table[i][j - 1] - affine_penalty - gap_penalty, upper_table[i][j - 1] - affine_penalty - gap_penalty)  # the value of moving from left to current position in dp table
                left_row.append(left)
                # for the E table
                up = max(upper_table[i - 1][j] - gap_penalty, mid_table[i - 1][j] - affine_penalty - gap_penalty, left_table[i - 1][j] - affine_penalty - gap_penalty)   # the value of moving from upper to current position in dp table
                upper_row.append(up)
                # for the D table
                seq = seq1[i - 1] + seq2[j - 1]
                score = subst_dict[seq]
                diagonal = mid_table[i - 1][j - 1] + score

                maximum = max(left, up, diagonal)
                mid_row.append(maximum)
                # can take more than one direction, so as long as it is equal to max value, append it to the list at position of pointer table
                if left == maximum:
                    pointer_row[j].append("left")
                if up == maximum:
                    pointer_row[j].append("up")
                if diagonal == maximum:
                    pointer_row[j].append("diagonal")
            # the first column
            elif i > 0 and j == 0:
                upper_row.append(-affine_penalty-gap_penalty*i)
                left_row.append(-10000000)
                mid_row.append(max(-10000000,-affine_penalty-gap_penalty*i))
                pointer_row[j].append("up")
            # the first row
            elif i == 0 and j > 0:
                upper_row.append(-10000000)
                left_row.append(-affine_penalty-gap_penalty*j)
                mid_row.append(max(-10000000,-affine_penalty-gap_penalty*j))
                pointer_row[j].append("left")
            # the [0,0] position, storing 0 for dp and "None" for pointer
            elif i == 0 and j == 0:
                mid_row.append(0)
                upper_row.append(-affine_penalty)
                left_row.append(-affine_penalty)
                pointer_row[j].append("None")
            mid_table[i] = mid_row
            upper_table[i] = upper_row
            left_table[i] = left_row
            pointer_table[i] = pointer_row
    nrow = len(seq1)
    ncol = len(seq2)
    bottom1 = ''
    bottom2 = ''
    top1 = ''
    top2 = ''
    # traverse back from the last position of table to find the paths and sequences of optimal alignments
    while nrow != 0 or ncol != 0:
        char1 = seq1[nrow - 1]
        char2 = seq2[ncol - 1]
        pointer = pointer_table[nrow][ncol]
        # top-most
        if "up" in pointer:
            nrow -= 1
            top1 += char1
            top2 += '-'
        elif "diagonal" in pointer:
            nrow -= 1
            ncol -= 1
            top1 += char1
            top2 += char2
        elif "left" in pointer:
            ncol -= 1
            top1 += '-'
            top2 += char2

    nrow = len(seq1)
    ncol = len(seq2)
    while nrow != 0 or ncol != 0:
        char1 = seq1[nrow - 1]
        char2 = seq2[ncol - 1]
        pointer = pointer_table[nrow][ncol]
        # bottom-most
        if "left" in pointer:
            ncol -= 1
            bottom1 += '-'
            bottom2 += char2
        elif "diagonal" in pointer:
            nrow -= 1
            ncol -= 1
            bottom1 += char1
            bottom2 += char2
        elif "up" in pointer:
            nrow -= 1
            bottom1 += char1
            bottom2 += '-'
    return (mid_table[len(seq1)][len(seq2)], (top1[::-1], top2[::-1]), (bottom1[::-1], bottom2[::-1]))


if __name__ == '__main__':
    run_ag_aligner_plus()
