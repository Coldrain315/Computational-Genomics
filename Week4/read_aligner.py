from collections import Counter
from bwt import *
from compsci260lib import *


def find(query, bwt_data):
    """Given a query sequence and a series of data structures containing
    various information about the reference genome, return a list containing
    all the locations of the query sequence in the reference genome.

    Args:
        query (str): query sequence to identify in the reference genome
        bwt_data (tuple): BWT data structures containing information about
        the reference genome. Refer to the documentation for `make_all` for
        more information.

    Returns:
        (list of ints): of all locations of the query sequence in the
                        reference genome
    """

    bwt, suffix_array, ranks, counts = bwt_data
    char_index = 0
    for char in query:
        # If there is no such element in the rank keys, I will directly return an empty set
        if char not in ranks.keys():
            return []

        rank_char = ranks[char]
        # get the number of elements and the location of each element
        count_elements = list(set(rank_char))
        element_index_lst = list()
        for count_element in count_elements:
            # skip when the count_element equals 0
            if count_element == 0:
                continue
            # get all the element locations in the suffix array of the reference genome
            element_index_lst.append(suffix_array[rank_char.index(count_element)] - 1 - char_index)
        # If it is the start position store it in the result list
        if char_index == 0:
            results = element_index_lst
        else:
            results = list(set(results) & set(element_index_lst))
        char_index += 1
    return sorted(results)


# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.

def rank(bwt_seq):
    """Takes as input a string transformed by the BWT. Returns a dictionary
    with characters as keys and lists as values. Each list contains the total
    number of occurrences for the corresponding character up until each
    position in the BWT-transformed string (i.e., its rank).

    For example:
        rank('ACTGA$TA')['$'] --> [0, 0, 0, 0, 0, 1, 1, 1]
        rank('ACTGA$TA')['A'] --> [1, 1, 1, 1, 2, 2, 2, 3]
        rank('ACTGA$TA')['C'] --> [0, 1, 1, 1, 1, 1, 1, 1]
        rank('ACTGA$TA')['G'] --> [0, 0, 0, 1, 1, 1, 1, 1]
        rank('ACTGA$TA')['T'] --> [0, 0, 1, 1, 1, 1, 2, 2]
    """
    rank = {}
    characters = set(bwt_seq)
    for character in characters:
        rank[character] = [0]
    rank[bwt_seq[0]] = [1]
    for letter in bwt_seq[1:]:
        for k, v in list(rank.items()):
            v.append(v[-1] + (k == letter))
    return rank


def make_suffix_array(seq):
    """Returns the suffix array of the input string.

    For example: print make_suffix_array('GATTACA$') --> [7, 6, 4, 1, 5, 0, 3, 2]
    """
    suffixes = {}
    for x in range(len(seq)):
        suffixes[seq[x:]] = x
    suffix_array = [suffixes[suffix] for suffix in sorted(suffixes.keys())]
    return suffix_array


def count_smaller_chars(seq):
    """Takes as input a string. Returns a dictionary with characters as keys
    and integers as values. The integers track the number of characters in the
    input string which are lexicographically smaller than the corresponding
    character key.

    For example, using an input DNA sequence like 'GATTACA':
        count_smaller_chars('GATTACA')['A'] --> 0
            (A, being lexicographically first in a DNA sequence,
            should always return 0)

        count_smaller_chars('GATTACA')['C'] --> 3
            (C, being second, should return the number of A's, which here is 3)

        count_smaller_chars('GATTACA')['G'] --> 4
            (G, being third, should return the number of A's or C's,
            which here is 4)

        count_smaller_chars('GATTACA')['T'] --> 5
            (T, being fourth, should return the number of A's or C's or G's,
            which here is 5)
    """
    characters = set(seq)
    cntr = Counter(seq)
    total = 0
    counts = {}
    for character in sorted(characters):
        counts[character] = total
        total += cntr[character]
    return counts


def make_all(reference):
    """Takes as input a reference string.

    Returns the data structures necessary to perform efficient exact
    string matching searches.
    """
    counts = count_smaller_chars(reference)
    reference = reference + '$'
    suffix_array = make_suffix_array(reference)
    bwt = forward_bwt(reference)
    ranks = rank(bwt)
    return bwt, suffix_array, ranks, counts


if __name__ == '__main__':
    ref = make_all("GATTACGATACTACG")
    # print(make_all("GATTACATACTACG"))
    print(find('TACG', ref))



