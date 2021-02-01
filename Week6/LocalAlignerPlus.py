
import sys, re
from compsci260lib import *
from aligner_helpers import *


def run_local_aligner_plus():
    """
    Locally align O18381.fasta and P63015.fasta and report the optimal score 
    and optimal local alignment(s) information, depending on whether you are 
    using Option A or B.
    """

    O18381 = get_fasta_dict('O18381.fasta')['sp|O18381|PAX6_DROME Paired box protein Pax-6 OS=Drosophila melanogaster OX=7227 GN=ey PE=1 SV=3']
    P63015 = get_fasta_dict('P63015.fasta')['sp|P63015|PAX6_MOUSE Paired box protein Pax-6 OS=Mus musculus OX=10090 GN=Pax6 PE=1 SV=1']
    # seq1 = "TAG"
    # seq2 = "TAAGATAAG"

    match = 2
    mismatch = -1
    gap_penalty = 8

    seq_type = validate_sequences(O18381, P63015)
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

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output_tuple = solve_local_aligner_plus(O18381, P63015, subst_dict, gap_penalty)
    print("The optimal local alignment score is %d" % output_tuple[0])
    print("The total number of locations in the table achieving this optimal score is %d" % len(output_tuple[1]))
    print("The locations in the table of each of these optima are", output_tuple[1])
    print("One of the optimal alignments is %s, from %d to %d on sequence1 and %s, from %d to %d on sequence2" %(output_tuple[2][-2],output_tuple[2][0][0],output_tuple[2][0][1],output_tuple[2][-1],output_tuple[2][1][0],output_tuple[2][1][1]))


def solve_local_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The procedure for collecting the inputs, running the aligner,
       and returning the score and optimal alignment(s).

       Note that for each returned local alignment, starting positions also 
       need to be returned. These are the positions of the first character in
       each aligned sequence relative to the original sequence.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character); this
        value should be positive because we will subtract it

    There are two possible implementations (and return values) for this
    problem (see the problem set instructions).

    Option A:

        A max score may be in multiple locations, so return the optimal score,
        the locations of all the maxima, and any one optimal alignment as a
        tuple.

        Returns a tuple of:
            (the optimal alignment score as an int, 

             the locations of the maxima in the dp table as a list of tuples.
             these positions will include the offset of the initialized penalty
             row and column, so that location (i,j) refers to the i-prefix of X
             and the j-prefix of Y, just as in lecture,

             tuple for an optimal alignment)

            The alignment will be in the form:

                  (tuple of indices of the characters of the first aligned 
                   sequence used in the alignment)

                  (tuple of indices of the characters of the second aligned 
                   sequence used in the alignment)

                  the first aligned sequence as a string,
                  the second aligned sequence as a string)

            As an example with the sequences:

                Sequence 1: TAG
                Sequence 2: TAAGATAAG

            A possible return may be:

                (11, # the optimal score

                 # the two maximal locations in the dp table
                 [(3, 4), (3, 9)],

                 # one possible alignment:
                 ((1, 3), # the nt positions mapping TA-G from TAG
                  (1, 4), # the nt positions mapping TAAG from TAAGATAAG

                  'TA-G', 'TAAG') # the sequences of the alignment
                )

            Corresponding to two maxima with an optimal
            alignment score of 11.

    Option B:

       Perform a top-most and bottom-most traceback for each of these max 
       scores. And return each in the form specified below.

        Returns a tuple of:
            (the optimal alignment score as an int, 

             tuple for top-most alignment from first traceback,
             tuple for bottom-most alignment from first traceback,

             ...
             
             tuple for top-most alignment from last traceback
             tuple for bottom-most alignment from last traceback
            )

            Each alignment will be in the form:
                 (the starting position of the first aligned sequence 
                    (0-indexed),
                  the starting position of the second aligned sequence
                    (0-indexed),
                  the first aligned sequence as a string,
                  the second aligned sequence as a string)

            As an example with the sequences:

                Sequence 1: TAG
                Sequence 2: TAAGATAAG

            A possible return may be:

                (11, 
                 (0, 0, 'T-AG', 'TAAG'), # first traceback
                 (0, 0, 'TA-G', 'TAAG'),
                 
                 (0, 5, 'T-AG', 'TAAG'), # second traceback
                 (0, 5, 'TA-G', 'TAAG'))

            Corresponding to two tracebacks each containing distinct top-most 
            and bottom-most alignments all scoring 11.
    """

    #
    # YOUR CODE GOES HERE
    #
    # Option A:
    # return (-1, [(-1, -1)], ((-1, -1), (-1, -1), '', ''))
    # Option B:
    # return (-1, (-1, -1, '', ''), (-1, -1, '', ''))
    list1 = []
    # initialize two two-D array with n rows. and n is the length of seq1 + 1
    dp_table = [list(list1) for i in range(len(seq1) + 1)]
    pointer_table = [list(list1) for i in range(len(seq1) + 1)]
    cur_max = 0
    position = []
    # go through all positions of two sequences. Plus one because there is an extra list and a column for '-'
    for i in range(len(seq1) + 1):
        dp_row = [list(list1) for i in range(len(seq2) + 1)]   # initialize an empty row for dp table
        pointer_row = [list(list1) for i in range(len(seq2) + 1)]  # initialize an array for pointer table

        for j in range(len(seq2) + 1):
            if i > 0 and j > 0:
                # for scores gained at each position, compare them with 0, if larger, takes the value, otherwise takes 0
                left = max(0, dp_table[i][j - 1] - gap_penalty)  # the value of moving from left to current position in dp table
                up = max(0, dp_table[i - 1][j] - gap_penalty)  # the value of moving from upper to current position in dp table
                seq = seq1[i - 1] + seq2[j - 1]
                score = subst_dict[seq]
                diagonal = max(dp_table[i - 1][j - 1] + score, 0)

                maximum = max(left, up, diagonal)
                dp_row[j] = maximum

                # updating the maximum score; if current best is larger than all best, update all best and refresh position
                if maximum > cur_max:
                    cur_max = maximum
                    position = [(i,j)]
                # if current best is equal to all best, add new position info to the list because we want locations of all best scores
                elif maximum == cur_max:
                    position.append((i,j))
                # can take more than one direction, so as long as it is equal to max value, append it to the list at position of pointer table
                if left == maximum:
                    pointer_row[j].append("left")
                if up == maximum:
                    pointer_row[j].append("up")
                if diagonal == maximum:
                    pointer_row[j].append("diagonal")
            # the first column, storing i*-gap for dp and "up" for pointer
            else:
                dp_row[j] = 0
                pointer_row[j] = "None"

            dp_table[i] = dp_row
            pointer_table[i] = pointer_row
    # take the first position with maximum score
    nrow = position[0][0]
    ncol = position[0][1]
    top1 = ''
    top2 = ''
    while dp_table[nrow][ncol] != 0:
        char1 = seq1[nrow - 1]
        char2 = seq2[ncol - 1]
        pointer = pointer_table[nrow][ncol]
        # bottom-most
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
    # display_dp_table(seq1, seq2, dp_table)
    return (cur_max, position, ((nrow+1, position[0][0]), (ncol+1, position[0][1]), top1[::-1], top2[::-1]))


if __name__ == '__main__':
    run_local_aligner_plus()
