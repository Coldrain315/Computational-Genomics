import sys
from compsci260lib import *


def run_global_aligner():
    """
    Given two sequences of either DNA or amino acids, initialize the
    appropriate substitution matrix, run the global aligner and report the
    optimal alignment score.
    """

    seq1 = "TAGTA"
    seq2 = "GAATCGGA"
    match = 2
    mismatch = -1
    gap_penalty = 2

    seq_type = validate_sequences(seq1, seq2)
    if seq_type == 1:
        # Both the sequences are DNA sequences so use the scores for match and
        # mismatch
        subst_matrix = create_subst_matrix_dna(match, mismatch)
        print(subst_matrix)
    elif seq_type == 2:
        # Both the sequences are protein sequences so read in the BLOSUM62
        # substitution matrix
        subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    else:
        sys.exit('Input sequences are of different types: not both DNA or both protein')

    # Obtain a dictionary of scores for aligning a pair of characters
    subst_dict = create_subst_matrix_dict(subst_matrix)
    optimal_score = solve_global_aligner(seq1, seq2, subst_dict, gap_penalty)
    print()
    print("The optimal score is %d." % optimal_score)
    # Report the optimal global alignment score
    # 
    # YOUR CODE GOES HERE
    #


def solve_global_aligner(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    and displaying the table and return the optimal alignment score

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): penalty for a column containing a gap char (g: use a 
            positive value because this value will be subtracted)

    Returns:
        (int) the optimal alignment score
    """

    # Initialize the DP table's data structure
    # as a list of lists of ints
    #
    # YOUR CODE HERE
    #
    list1 = []
    # initialize a two-D array with i rows. and i is the length of seq1
    dp_table = [list(list1) for i in range(len(seq1) + 1)]
    # go through all positions of two sequences. Plus one because there is an extra list and a column for '-'
    for i in range(len(seq1) + 1):
        dp_row = []  # initialize an empty row for dp table
        for j in range(len(seq2) + 1):
            if i > 0 and j > 0:
                left = dp_table[i][j - 1] - gap_penalty  # the value of moving from left to current position in dp table
                up = dp_table[i - 1][j] - gap_penalty  # the value of moving from upper to current position in dp table

                seq = seq1[i - 1] + seq2[j - 1]  # combining
                score = subst_dict[seq]  # check the dictionary that contains score for combinations of chars
                diagonal = dp_table[i - 1][j - 1] + score
                maximum = max(left, up, diagonal)
                dp_row.append(maximum)
            elif i > 0 and j == 0:
                dp_row.append(-gap_penalty * i)
            elif i == 0 and j > 0:
                dp_row.append(-gap_penalty * j)
            elif i == 0 and j == 0:
                dp_row.append(0)
            dp_table[i] = dp_row

    # Compute the score of the optimal global alignment
    max_value = global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table)

    # Display the dp table
    display_dp_table(seq1, seq2, dp_table)
    return max_value


def create_subst_matrix_dna(match, mismatch):
    """Create a substitution matrix for DNA sequences. This method
    should be used with create_subst_matrix_dict().

    Args:
        match (int): score for a match
        mismatch (int): score for a mismatch

    Returns:
        (list of strings): representing the substitution matrix for DNA. Scores
        are separated by a single space.
    """
    subst_matrix = []
    subst_matrix.append("A C G T")
    subst_matrix.append("%s %s %s %s" % (match, mismatch, mismatch, mismatch))
    subst_matrix.append("%s %s %s %s" % (mismatch, match, mismatch, mismatch))
    subst_matrix.append("%s %s %s %s" % (mismatch, mismatch, match, mismatch))
    subst_matrix.append("%s %s %s %s" % (mismatch, mismatch, mismatch, match))

    return subst_matrix


def create_subst_matrix_aa(blosum_path):
    """Create a substitution matrix for amino acid sequences. This method
    should be used with create_subst_matrix_dict().

        Args:
            blosum_path (str): path to the BLOSUM file to load

        Returns:
            (list of strings): representing the substitution matrix for amino
            acids. Scores in each row are separated by spaces.
        """

    subst_matrix = load_matrix(blosum_path)
    return subst_matrix


def global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table):
    """A dynamic programming algorithm that takes two sequences and returns the
    score of the optimal alignment.

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score
    """

    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)      # so I is 1 more than m
    J = len(dp_table[0])   # so J is 1 more than n

    # Initialize the dp table with solutions to base cases using linear gap
    # penalty
    #
    # YOUR CODE GOES HERE
    #
    for cols in range(J):
        dp_table[0][cols] = -cols*gap_penalty
    for rows in range(I):
        dp_table[rows][0] = -rows*gap_penalty

    # Compute the scores for the rest of the matrix,
    # i.e. all the elements in dp_table[i][j] for i > 0 and j > 0.
    #
    # YOUR CODE GOES HERE
    #
    for rows in range(1, I):
        for cols in range(1, J):
            dp_table[rows][cols] = max([dp_table[rows-1][cols-1]+subst_dict[seq1[rows-1]+seq2[cols-1]],dp_table[rows-1][cols]-gap_penalty,dp_table[rows][cols-1]-gap_penalty])
    # The optimal score is found at the lower right corner of the dp table:
    return dp_table[I-1][J-1]


def compute_global_aligner_score_linspace(seq1, seq2, subst_dict, gap_penalty):
    """Extra Challenge: compute the score (not the actual alignment) of the
       best global alignment in O(mn) time but only O(min(m,n)) space

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score
    """

    seqshort, seqlong = (seq1, seq2) if len(seq1) < len(seq2) else (seq2, seq1)

    #
    # YOUR CODE GOES HERE
    #
    return -1


########################################################################
########################################################################
# Below are some complete functions that may be handy to call
########################################################################
########################################################################


def display_dp_table(seq1, seq2, dp_table):
    """A function that displays the dp table once it's been computed."""
    # Add gap character in front of the two sequences
    seq1 = "-" + seq1
    seq2 = "-" + seq2

    # Trying to get the columns lined up; this works as long as |values| < 99,
    # and can be modified for bigger values
    print("\n\n\t\t ", end=' ')

    for i in range(len(seq2)):
        print(seq2[i] + "\t ", end=' ')

    for i in range(len(seq1)):
        print("\n\t" + seq1[i], end=' ')

        for j in range(len(seq2)):
            adjustment = '  '
            if dp_table[i][j] > 9:
                adjustment = ' '
            if dp_table[i][j] < 0:
                adjustment = ' '
            if dp_table[i][j] < -9:
                adjustment = ''
            print("\t", adjustment, dp_table[i][j], end=' ')


def is_dna(in_str):
    """Determine if string is a valid DNA sequence; it's valid if it only
    contains {ACGT} characters and is not empty."""
    if re.findall(r"^[ACGT]+$", in_str.upper()):
        return True
    else:
        return False


def is_protein(in_str):
    """Determine if string is a valid protein sequence; it's valid if it only
    contains {ACGTRNDQEHILKMFPSWYV} characters and is not empty."""
    if re.findall(r"^[ACGTRNDQEHILKMFPSWYV]+$", in_str):
        return True
    else:
        return False


def validate_sequences(seq1, seq2):
    """Determines if a given pair of sequences are both valid DNA sequences,
    both valid protein sequences, or an invalid pair of sequences."""
    if is_dna(seq1) and is_dna(seq2):
        return 1
    elif is_protein(seq1) and is_protein(seq2):
        return 2
    else:
        return 0


def get_sequence():
    """A console interface for sequence input.

    Returns a string sequence the user entered or from the file the user
    specified.
    """
    print("\nChoose one of the following options")
    print("\n\t1 The sequence will be entered by hand")
    print("\n\t2 The sequence will be read from a file")
    print("\nEnter choice: ")

    while True:
        try:
            choice = int(input())
            break
        except Exception:
            print("Please enter an integer: ")

    seq = ""
    if choice == 1:
        print("\nEnter a sequence: ")
        sys.stdout.flush()
        seq = input()
        seq = seq.strip().upper()
    else:
        print("\nEnter the filename: ")
        sys.stdout.flush()
        file_name = input()
        file_name = file_name.strip()
        seq = list(get_fasta_dict(file_name).values())[0]

    return seq


def create_subst_matrix_dict(rows):
    """Reads in a list of rows of strings.

    The first row labels the character types (nucleotide types if DNA,
    amino acid types if protein). The remaining rows are strings that
    contain the values in the substitution matrix.  Function transforms
    these rows of strings into a matrix, and then transforms that into a
    dictionary where the key is a pair of characters (literally, a
    string of length two) and the value is the score of aligning those
    characters. The reason for doing this is so that we can index into a
    data structure using the characters themselves.  Function then
    returns the dictionary.  Gap score is not included in any of this;
    that's handled separately.

    Args:
        (list of strings): representation of the substitution matrix

    Returns:
        (dictionary string -> int): dictionary representation of the
        substitution matrix mapping a string representing the row and column
        characters to an alignment score of the two characters:

            <first character> + <second character> -> <alignment score>

        The returned dictionary may look something like this:
        {
            'AA': 1,
            'AT': 0,
            ...
        }

    """
    matrix = []

    # Remove whitespace in front of the first row
    rows[0].strip()

    # The first row of the matrix contains a sequence of characters.  Store
    # them as a list. These characters are labels for the rows & columns of the
    # substitution matrix.
    reference_list = rows[0].split()

    # Now get rid of the first row
    rows.pop(0)

    for row in rows:
        # Again, remove whitespace in front of rows
        row = row.strip()
        cols = row.split()
        matrix.append(cols)

    # Create the dictionary
    subst_dict = {}
    for i in range(len(reference_list)):
        for j in range(len(reference_list)):
            # In the following line, the + is used to join the two characters
            # into a string of length 2.
            subst_dict[reference_list[i] +
                       reference_list[j]] = int(matrix[i][j])

    return subst_dict


def load_matrix(in_file):
    """Function to read a substitution matrix in from a file."""

    with open(in_file, 'r') as f:
        # Check to see if the file exists
        while f is None:
            print("File does not exist. \nEnter the file name again: ")
            in_file = input()
            in_file = in_file.strip()
            f = open(in_file, 'r')

        lines = f.readlines()

    return lines


if __name__ == '__main__':
    run_global_aligner()
