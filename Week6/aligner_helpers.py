
import re

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


def create_subst_matrix_dna(match,mismatch):
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

