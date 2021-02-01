import sys, re
# import numpy as np
from compsci260lib import *
from GlobalAligner import *


def run_global_aligner_plus():
    """Generate the optimal global alignments between:

        2.d)
        atpa_Hs.fasta, atpaEc.fasta

        2.f)
        atpaMm.fasta, atpaHs.fasta
        atpaMm.fasta, atpaEc.fasta 

        atpaBs.fasta, atpaHs.fasta
        atpaBs.fasta, atpaEc.fasta 

        atpaMm.fasta, atpaBs.fasta

    For each alignment, report the optimal alignment score, 
    the top-most alignment, and the bottom-most alignment.
    """

    #
    # YOUR CODE GOES HERE
    # Bs = get_fasta_dict('atpa_Bs.fasta')[
    #     'sp|P37808|ATPA_BACSU ATP synthase subunit alpha OS=Bacillus subtilis (strain 168) OX=224308 GN=atpA PE=1 SV=3']
    # Mm = get_fasta_dict('atpa_Mm.fasta')[
    #     'sp|Q03265|ATPA_MOUSE ATP synthase subunit alpha, mitochondrial OS=Mus musculus OX=10090 GN=Atp5f1a PE=1 SV=1']
    Ec = get_fasta_dict('atpa_Ec.fasta')[
        'sp|P0ABB0|ATPA_ECOLI ATP synthase subunit alpha OS=Escherichia coli (strain K12) OX=83333 GN=atpA PE=1 SV=1']
    Hs = get_fasta_dict('atpa_Hs.fasta')[
        'sp|P25705|ATPA_HUMAN ATP synthase subunit alpha, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1A PE=1 SV=1']
    # print(protein_try2)
    seq1 = "TAGTA"
    seq2 = "GAATCGGA"
    # seq1 = "ATAGG"
    # seq2 = "ATCCGG"
    match = 2
    mismatch = -1
    gap_penalty = 8
    # gap_penalty = 2

    seq_type = validate_sequences(Ec, Hs)
    # seq_type = validate_sequences(seq1, seq2)

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
    # solve_global_aligner_plus(Ec, Hs, subst_dict, gap_penalty)

    # c)
    output_tuple = solve_global_aligner_plus(Hs, Ec, subst_dict, gap_penalty)
    print("The optimal alignment score for comparison between atpa_Hs and atpa_Ec is %d." % output_tuple[0])
    print("The top most alignment achieving this score is:")
    writing_aligners_helper(output_tuple[1][0], output_tuple[1][1], char_line=50)
    print("The bottom most alignment achieving this score is:")
    writing_aligners_helper(output_tuple[2][0], output_tuple[2][1], char_line=50)
    #
    # # f)
    # list_of_sequence = {"Bs": Bs, "Mm": Mm, "Ec": Ec, "Hs": Hs}
    # key_list = list(list_of_sequence.keys())
    # for i in range(len(key_list)):
    #     for j in range(i+1,len(key_list)):
    #         if key_list[i] == "Ec" and key_list[j] == "Hs":
    #             continue
    #         else:
    #             subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    #             subst_dict = create_subst_matrix_dict(subst_matrix)
    #             output_tuple = solve_global_aligner_plus(list_of_sequence[key_list[i]], list_of_sequence[key_list[j]], subst_dict, gap_penalty)
    #             print("The optimal alignment score for comparison between atpa_%s and atpa_%s is %d." % (key_list[i],key_list[j],output_tuple[0]))
    #             print("The top most alignment achieving this score is:")
    #             writing_aligners_helper(output_tuple[1][0], output_tuple[1][1], char_line=50)
    #


def solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    and displaying the table and final value.

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

    Returns a tuple of:
        (the optimal alignment score as an int, 
         the top-most alignment achieving this score as a tuple of strings,
         the bottom-most alignment achieving this score as a tuple of strings)

        Example output:

            (6, ("AT-AGG", "ATCCGG"), ("ATA-GG", "ATCCGG"))

    Note: If you do the extra challenge to report all optimal alignments, 
    you can lengthen the size of the return tuple, but ensure that its second
    element (i.e., the first alignment) remains the top-most alignment, while 
    the last element is the bottom-most alignment. e.g.

        (optimal score, (top-most alignment sequences), ..., 
         (bottom-most alignment sequences))

    """

    #
    list1 = []
    # initialize two two-D array with n rows. and n is the length of seq1 + 1
    dp_table = [list(list1) for i in range(len(seq1) + 1)]
    pointer_table = [list(list1) for i in range(len(seq1) + 1)]
    # go through all positions of two sequences. Plus one because there is an extra list and a column for '-'
    for i in range(len(seq1) + 1):
        dp_row = []  # initialize an empty row for dp table
        pointer_row = [list(list1) for i in range(len(seq2) + 1)]  # initialize an array for pointer table
        for j in range(len(seq2) + 1):
            if i > 0 and j > 0:
                left = dp_table[i][j - 1] - gap_penalty  # the value of moving from left to current position in dp table
                up = dp_table[i - 1][j] - gap_penalty  # the value of moving from upper to current position in dp table

                seq = seq1[i - 1] + seq2[j - 1]
                score = subst_dict[seq]
                diagonal = dp_table[i - 1][j - 1] + score
                maximum = max(left, up, diagonal)
                dp_row.append(maximum)
                # can take more than one direction, so as long as it is equal to max value, append it to the list at position of pointer table
                if left == maximum:
                    pointer_row[j].append("left")
                if up == maximum:
                    pointer_row[j].append("up")
                if diagonal == maximum:
                    pointer_row[j].append("diagonal")
            # the first column, storing i*-gap for dp and "up" for pointer
            elif i > 0 and j == 0:
                dp_row.append(-gap_penalty * i)
                pointer_row[j].append("up")
            # the first row, storing j*-gap for dp and "left" for pointer
            elif i == 0 and j > 0:
                dp_row.append(-gap_penalty * j)
                pointer_row[j].append("left")
            # the [0,0] position, storing 0 for dp and "None" for pointer
            elif i == 0 and j == 0:
                dp_row.append(0)
                pointer_row[j].append("None")
            dp_table[i] = dp_row
            pointer_table[i] = pointer_row
    max_value = global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table)
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
    return (max_value, (top1[::-1], top2[::-1]), (bottom1[::-1], bottom2[::-1]))


def writing_aligners_helper(seq1, seq2, char_line=50):
    """
    The help function for separating the sequence in different chunks.
    Args:
        seq1(str): The aligned sequence 1 with "-"s
        seq2(str): The aligned sequence 2 with "-"s
        char_line(int): The number of character to print each line
    """
    # get the length of the sequence
    len_seq = len(seq1)
    # calculate the number of line it need if each line have (char_line) characters
    len_line = len_seq // char_line + 1

    for index in range(len_line):
        print(seq1[index * char_line:(index + 1) * char_line])
        print(seq2[index * char_line:(index + 1) * char_line])
        print()


if __name__ == '__main__':
    run_global_aligner_plus()
