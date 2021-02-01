from compsci260lib import *
from UltrametricAdditive import is_ultrametric, is_additive
from GlobalAlignerPlus import *


def build_tree():
    """
    Read the aligned student sequences, compute the dictionary of distances,
    construct the appropriate phylogenetic tree.
    """

    # Load the aligned student sequences and create the dictionary of distances
    # in the same format as found in UltrametricAdditive.py. You will need to
    # convert the student names to integers ("Student_1" -> 1)
    students_alignments = get_fasta_dict("students.aligned.fasta")
    # print(students_alignments.keys())
    # YOUR CODE HERE
    new_sa = students_alignments
    for i in list(students_alignments.keys()):
        new_key = i.split('_')[1]
        new_sa[new_key] = new_sa.pop(i)
    dist = {}
    for key1 in range(len(list(new_sa.keys()))):
        for key2 in range(key1 + 1, len(list(new_sa.keys()))):
            key = str(key1 + 1) + ',' + str(key2 + 1)
            value = compute_dist(new_sa[str(key1 + 1)], new_sa[str(key2 + 1)])
            dist[key] = value
    print("The table of pairwise distances is:")
    print_dist(dist)
    dist_1 = {"1,2": 0.8, "1,3": 0.4, "1,4": 0.6, "1,5": 0.8,
              "2,3": 0.8, "2,4": 0.8, "2,5": 0.4,
              "3,4": 0.6, "3,5": 0.8,
              "4,5": 0.8}

    # Check if the distances are ultrametric
    threshold = 0.003  # problem-specific distance threshold for this problem
    if is_ultrametric(dist, threshold=threshold):
        print("\nThe distance is ultrametric.")
    else:
        print("\nThe distance is not ultrametric.")

    # Check if the distances are additive
    students_alignments = get_fasta_dict("students.aligned.fasta")
    new_sa = students_alignments
    for i in list(students_alignments.keys()):
        new_key = i.split('_')[1]
        new_sa[new_key] = new_sa.pop(i)
    dist = {}
    for key1 in range(len(list(new_sa.keys()))):
        for key2 in range(key1 + 1, len(list(new_sa.keys()))):
            key = str(key1 + 1) + ',' + str(key2 + 1)
            value = compute_dist(new_sa[str(key1 + 1)], new_sa[str(key2 + 1)])
            dist[key] = value
    if is_additive(dist, threshold=threshold):
        print("\nThe distance is additive.\n")
    else:
        print("\nThe distance is not additive.\n")

    # Report the Newick representation of the phylogenetic tree    
    #
    students_alignments = get_fasta_dict("students.aligned.fasta")
    new_sa = students_alignments
    for i in list(students_alignments.keys()):
        new_key = i.split('_')[1]
        new_sa[new_key] = new_sa.pop(i)
    dist = {}
    for key1 in range(len(list(new_sa.keys()))):
        for key2 in range(key1 + 1, len(list(new_sa.keys()))):
            key = str(key1 + 1) + ',' + str(key2 + 1)
            value = compute_dist(new_sa[str(key1 + 1)], new_sa[str(key2 + 1)])
            dist[key] = value
    result = compute_nj_newick([2,4,6,3,8,7,5,1], dist)
    print("The Newick format of reconstructed tree is", result)
    #

    # Call `summarize_alignment` to report the differences between the two
    # most similar and two most different sequence alignments
    #
    sars = [2, 3, 4, 6, 8]
    min_all = 0
    min_all_list = []
    min_cov = 0
    min_cov_list = []
    max_cov = 0
    max_cov_list = []
    for i in range(1, len(new_sa)+1):
        for j in range(i+1, len(new_sa) + 1):
            seq1 = new_sa[str(i)]
            seq2 = new_sa[str(j)]
            result = summarize_alignment(seq1, seq2)
            if result[1] > min_all:
                min_all = result[1]
                min_all_list = [(i, j)]
            elif result[1] == min_all:
                min_all_list.append((i, j))
            if i in sars and j in sars:
                if result[1] > min_cov:
                    min_cov = result[1]
                    min_cov_list = [(i, j)]
                elif result[1] == min_cov:
                    min_cov_list.append((i, j))
                if result[0] > max_cov:
                    max_cov = result[0]
                    max_cov_list = [(i, j)]
                elif result[0] == max_cov:
                    max_cov_list.append((i, j))
                    # the differences between the two
                    # SARS-CoV-2 genome sequences that are most similar, the two SARS-CoV-2 genome sequences that are most
                    # different, and the two genome sequences that are most different overall
    print("The highest score of match in SARS-CoV-2 genome sequences is %d, and the two SARS-CoV-2 genome sequences that are most similar are student %d and student %d." % (max_cov,max_cov_list[0][0], max_cov_list[0][1]))
    print("The highest score of mismatch in SARS-CoV-2 genome sequences is %d, and the two SARS-CoV-2 genome sequences that are most similar are student %d and student %d, student %d and student %d." % (min_cov,min_cov_list[0][0], min_cov_list[0][1],min_cov_list[1][0], min_cov_list[1][1]))
    print("The highest score of mismatch overall is %d, and the two genome sequences that are most different are student %d and student %d." % (
        min_all, min_all_list[0][0], min_all_list[0][1]))

    #


def compute_nj_newick(seq_names, dist):
    """
    Performs neighbor joining to construct the phylogenetic tree from a
    distance table.

    Args:
        seq_names (list of ints): representing the sequence names.
                e.g. [1, 2, ..] for ['Student_1', 'Student_2', ...]

        dist (dict of str to float): distance table mapping pairs of students 
            to float distances. Refer to UltrametricAdditive.py for examples of
            distance tables.

    Returns:
        the Newick representation of the phylogenetic tree as a string
    """

    # ------------- Implementation of the Neighbor Joining algorithm ----------
    # Keeping track of variable names:

    # dist              - dictionary containing the computed pair-wise
    #                     distances between nodes
    # node_list         - array containing the list of current nodes (L), which
    #                     gradually decreases in size while iterating to build 
    #                     the tree
    # avg_dist          - dictionary containing the averaged distances 
    #                     (r values) for all of the current nodes (those in L)
    # D                 - dictionary containing the 'adjusted' pairwise
    #                     distances between nodes
    # newick            - dictionary to maintain the Newick representation of
    #                     each leaf node or subtree

    # ------------- Initialization: -------------------------------------------

    node_list = seq_names  # L = the list of current nodes, ['2', '4', '6', '3', '8', '7', '5', '1']
    newick = {}
    for i in range(1, len(node_list) + 1):
        newick[i] = "" + str(i)
    avg_dist = {}  # averaged distances (r values) for all current nodes

    for i in range(1, len(node_list) + 1):

        avg_dist[i] = 0
        for j in range(1, len(node_list) + 1):

            if i != j:
                avg_dist[i] += dist["%d,%d" % (i, j)] if i < j else \
                    dist["%d,%d" % (j, i)]

        avg_dist[i] = avg_dist[i] / (len(node_list) - 2)

    max_node = len(node_list)  # the maximum key used for a node
    # print(max_node, avg_dist)
    # -------------- Iteration to build the tree --------------------

    # As long as the list contains more than 2 nodes, iterate 
    while len(node_list) > 2:

        # ---------- Begin your code -------------

        # Compute the 'adjusted' distances between nodes using the original
        # distances (from dist)
        # and averaged distances (from avg_dist)

        # Let D be the dict to contain 'adjusted' distances
        D = {}

        # Loop through each pair of nodes as entered in the dist dict
        # Use the entries from the avg_dist and calculate entries for the D
        # dict
        # print(dist)
        for i in dist.keys():
            first = i.split(',')[0]
            second = i.split(',')[1]
            new_value = dist[i] - avg_dist[int(first)] - avg_dist[int(second)]
            D[i] = new_value
        # Pick the pair i,j in node_list for which adjusted distance D_ij is
        # minimal.
        minimal = 0
        I = ''
        J = ''
        for key, value in D.items():
            if value < minimal:
                minimal = value
                I = key.split(',')[0]
                J = key.split(',')[1]
        # You may find the function two_min_in_dict helpful.
        (i, j) = int(I), int(J)  # Replace with your pair
        # Define a new node k and set dist[m,k] for all nodes m in node_list
        # as (dist[i,m] + dist[j,m] - dist[i,j]) / 2
        k = max_node + 1
        # print(k)
        # node_list.append(k)
        list_key = []
        for n in dist.keys():
            key1 = n.split(',')[0]
            key2 = n.split(',')[1]
            if key1 not in list_key:
                list_key.append(key1)
            if key2 not in list_key:
                list_key.append(key2)
        for m in list_key:
            if int(m) != i and int(m) != j:
                dist_ij = dist[str(i) + ',' + str(j)]
                if int(m) < k:
                    new_key = m + ',' + str(k)
                elif int(m) > k:
                    new_key = str(k) + ',' + m
                if i > int(m):
                    dist_im = dist[m + ',' + str(i)]
                elif int(m) > i:
                    dist_im = dist[str(i) + ',' + m]
                if j > int(m):
                    dist_jm = dist[m + ',' + str(j)]
                elif int(m) > j:
                    dist_jm = dist[str(j) + ',' + m]
                new_value = (dist_im + dist_jm - dist_ij) / 2
                dist[new_key] = new_value
                # dist[new_key] = new_value
        # max_node had been earlier set to the largest key used for a node
        k = max_node + 1
        max_node += 1
        # m = 0
        # for ind in range(len(node_list)):
        #     m = node_list[ind]
        # ---------- End your code -------------

        # Add the new node k to the Newick format representation of the tree
        # with edges of lengths 
        # dik = (dist[i,j] + avg_dist[i] - avg_dist[j])/2
        # djk = dist[i,j]-d[i,k], 
        # joining k to i and j

        d_ik = (dist["%d,%d" % (i, j)] + avg_dist[i] - avg_dist[j]) / 2
        d_jk = dist["%d,%d" % (i, j)] - d_ik
        newick[k] = "(" + newick[i] + ":" + "%.7f" % (d_ik) + "," \
                    + newick[j] + ":" + "%.7f" % (d_jk) + ")"

        # Remove i and j from node_list and add k
        temp = []
        for idx in range(0, len(node_list)):

            if node_list[idx] != i and node_list[idx] != j:
                temp.append(node_list[idx])

        temp.append(k)

        node_list = list(temp)

        # Update the r terms
        if len(node_list) > 2:
            avg_dist[k] = 0

            for ind in range(0, len(node_list) - 1):
                m = node_list[ind]
                # print(m)
                avg_dist[m] = avg_dist[m] * (len(node_list) - 1)
                avg_dist[m] -= dist["%d,%d" % (m, i)] if m < i \
                    else dist["%d,%d" % (i, m)]
                avg_dist[m] -= dist["%d,%d" % (m, j)] if m < j \
                    else dist["%d,%d" % (j, m)]
                avg_dist[m] += dist["%d,%d" % (m, k)]

                avg_dist[m] /= (len(node_list) - 2)
                avg_dist[k] += dist["%d,%d" % (m, k)]

            avg_dist[k] = avg_dist[k] / (len(node_list) - 2)

        # Remove any elements from the dict that contain nodes i or j
        delete_from_dict(dist, i, j)
        delete_from_dict(D, i, j)
        delete_from_dict(newick, i, j)
        delete_from_dict(avg_dist, i, j)

    # Return the Newick representation
    return ("(" + newick[node_list[0]] + ":" +
            "%.7f" % (list(dist.values())[0]) +
            "," + newick[node_list[1]] + ":0);\n")


def summarize_alignment(seq1, seq2):
    """
    Summarize the alignment between two sequences by computing the number of
    matches, mismatches, and gaps. This code will contain similar logic to 
    the provided `compute_dist` function.

    Note: that we performed multiple sequence alignment to obtain the 
    aligned sequences in students.aligned.fasta. So, for any pair of sequences, 
    you may find a gap at the same place in the two aligned sequences. Gaps 
    should only be counted if they are matched with a non-gap character 
    (ignore a gap aligned to a gap).

    Args:
        seq1 (str): the first sequence, extracted from a multiple sequence
                    alignment
        seq2 (str): the second sequence, extracted from a multiple sequence
                    alignment

    Returns:
        a tuple of the number of (matches, mismatches, gaps) between the two
        sequences
    """

    #
    # YOU CODE HERE
    #
    mismatch = 0
    match = 0
    gap = 0
    for i in range(0, len(seq1)):
        if seq1[i] == "-":
            gap += 1
        elif seq2[i] == "-":
            gap += 1
        elif seq1[i] == seq2[i]:
            match += 1
        else:
            mismatch += 1

    # return the number of matches, mismatches and gaps as a tuple
    return (match, mismatch, gap)


########################################################################
# Provided functions for this problem
########################################################################

def compute_dist(seq1, seq2):
    """Returns the distance between two sequences. The distance is computed
    as the ratio between the number of mismatches and the total number
    of matches or mismatches.

    Args:
        seq1 (string) - first sequence
        seq2 (string) - second sequence

    Returns:
        the ratio of mismatches over the total number of matches or mismatches
        as a float
    """

    mismatch = 0
    match_or_mismatch = 0

    for i in range(0, len(seq1)):
        if seq1[i] == "-":
            continue
        elif seq2[i] == "-":
            continue
        elif seq1[i] == seq2[i]:
            match_or_mismatch += 1
        else:
            mismatch += 1
            match_or_mismatch += 1

    return float(mismatch) / match_or_mismatch


def print_dist(dist):
    """
    Print the distance table
    """

    idx = [int(i.split(',')[0]) for i in list(dist.keys())]
    idx.extend([int(i.split(',')[1]) for i in list(dist.keys())])
    max_idx = max(idx)

    print("\t", end=' ')
    for col in range(2, max_idx + 1):
        print("{:>7}".format(col), end=' ')
    print()

    for row in range(1, max_idx):
        print('%d\t' % row, end=' ')
        for col in range(2, row + 1):
            print('       ', end=' ')
        for col in range(row + 1, max_idx + 1):
            print('{:>7.4f}'.format(dist[str(row) + ',' + str(col)]), end=' ')
        print()

    ########################################################################


# Functions used for Neighbor Joining Algorithm
########################################################################

def delete_from_dict(dictionary, i, j):
    """Deletes the dict entries with keys that contain i or j."""

    for k in list(dictionary.keys()):
        ks = [int(_) for _ in str(k).split(',')]
        if i in ks or j in ks:
            del dictionary[k]


def min_in_dict(wiki):
    """Returns the key associated with the minimum value in the dict."""
    import operator
    return min(iter(wiki.items()), key=operator.itemgetter(1))[0]


def two_min_in_dict(dictionary):
    import operator
    sorted_dict = sorted(iter(dictionary.items()), key=operator.itemgetter(1))
    element = sorted_dict[0][0]  # get the first element of the tuple
    (i, j) = element.split(",")
    return (int(i), int(j))


if __name__ == '__main__':
    build_tree()
