from compsci260lib import *

# Note that the '$' character will be used to designate the end of a given
# string.


def forward_bwt(seq):
    """forward_bwt(seq) takes as input a string containing the EOF character to
    which the BWT must be applied. The method should then return the result of
    the BWT on the input string.

    For example:
        forward_bwt('GATTACA$') --> 'ACTGA$TA'

    Args:
        seq (str): input string with an EOF character

    Returns:
        (str): the transformed string
    """

    # list for rows with seq order rotated
    list_of_transform = []
    for i in range(len(seq)):
        # for every new row, move the char at first position to the end
        new_seq = seq[i:] + seq[0:i]
        list_of_transform.append(new_seq)
    # sort the rows alphabetically
    new_list = sorted(list_of_transform)
    transformed = ""
    # get the last column. and it will be bwt
    for j in new_list:
        transformed += j[-1]
    return transformed


def reverse_bwt(seq):
    """reverse_bwt(seq) takes as input a string containing the EOF character to
    which the reverse of the BWT must be applied. The method should then return
    the result of the reversal on the input string.

    For example:
        reverse_bwt('ACTGA$TA') --> 'GATTACA$'

    Args:
        seq (str): input string with an EOF character

    Returns:
        (str): the transformed string
    """

    seq_list = list(seq)  # turn the bwt sequence into list; this list will store the rows of matrix
    new_list = list(seq)
    for i in range(len(seq)-1):
        sort_seq = sorted(new_list)
        for j in range(len(seq)):
            new_list[j] = seq_list[j] + sort_seq[j]
    # for row in the matrix, and the row ends with $ is the original sequence
    for a in sorted(new_list):
        if a[-1] == "$":
            original = a

    return original


def solve_bwt():
    """Load the mystery.txt file and decode its contents using the 
    reverse BWT transformation. Report the decoded contents."""

    # Example sequences for forward_bwt() and reverse_bwt();
    # you may change them to different sequences to test your code.
    """
    seq1 = 'GATTACA$'
    seq2 = 'ACTGA$TA'
    forward_bwt(seq1)
    reverse_bwt(seq2)
    """

    # Code to open the file mystery.txt in the correct encoding
    # across platforms, and read its contents into a variable.
    """
    with open('mystery.txt', 'r', encoding='UTF-8') as f:
         mystery_seq = f.read()
    """

    seq = 'CGGACTAACGGACTAACGGACTAACGGACTAA$'
    bwt_seq = forward_bwt(seq)
    print("The result of burrow's wheeler transformed sequence is %s." % bwt_seq)
    file = open('mystery.txt', encoding='utf-8')
    text = file.read()
    reverse_bwt_text = reverse_bwt(text)
    reverse_bwt_text = reverse_bwt_text.replace('_',' ')
    print(reverse_bwt_text)


if __name__ == '__main__':
    solve_bwt()
