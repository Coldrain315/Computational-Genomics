from compsci260lib import *


def solve_ultrametric_additive():
    # Distance metrics for table 1 and table 2
    dist_1 = {"1,2": 0.3, "1,3": 0.7, "1,4": 0.9,
              "2,3": 0.6, "2,4": 0.8,
              "3,4": 0.6}  # fill in table 1 here

    dist_2 = {"1,2": 0.8, "1,3": 0.4, "1,4": 0.6, "1,5": 0.8,
              "2,3": 0.8, "2,4": 0.8, "2,5": 0.4,
              "3,4": 0.6, "3,5": 0.8,
              "4,5": 0.8}

    # Check if dist_1 and dist_2 are ultrametric and additive by
    # calling is_ultrametric and is_additive with the default
    # threshold value (1e-4).
    #
    result_ult1 = is_ultrametric(dist_1)
    print("The ultrametricity of D1 is", result_ult1)
    dist_1 = {"1,2": 0.3, "1,3": 0.7, "1,4": 0.9,
              "2,3": 0.6, "2,4": 0.8,
              "3,4": 0.6}  # fill in table 1 here
    result_add1 = is_additive(dist_1)
    print("The additivity of D1 is", result_add1)
    dict_ult2 = dist_2
    result_ult2 = is_ultrametric(dist_2)
    print("The ultrametricity of D2 is", result_ult2)
    dist_2 = {"1,2": 0.8, "1,3": 0.4, "1,4": 0.6, "1,5": 0.8,
              "2,3": 0.8, "2,4": 0.8, "2,5": 0.4,
              "3,4": 0.6, "3,5": 0.8,
              "4,5": 0.8}
    result_add2 = is_additive(dist_2)
    print("The additivity of D2 is", result_add2)

    # Construct the ATPA synthase distance metric table
    atpa_table = {"1,2":0.5,"1,3":0.5,"1,4":0.1,"1,5":0.4,"1,6":0.4,
                  "2,3":0.3,"2,4":0.5,"2,5":0.5,"2,6":0.5,
                  "3,4":0.5,"3,5":0.5,"3,6":0.5,
                  "4,5":0.4,"4,6":0.4,
                  "5,6":0.3}  # fill in ATPA synthase distance metric table

    # Determine if the ATPA synthase distance metrics
    # are ultrametric and additive using the default
    # threshold value (1e-4).
    #
    result_ult3 = is_ultrametric(atpa_table)
    print("The ultrametricity of ATPA synthase is", result_ult3)
    atpa_table = {"1,2": 0.5, "1,3": 0.5, "1,4": 0.1, "1,5": 0.4, "1,6": 0.4,
                  "2,3": 0.3, "2,4": 0.5, "2,5": 0.5, "2,6": 0.5,
                  "3,4": 0.5, "3,5": 0.5, "3,6": 0.5,
                  "4,5": 0.4, "4,6": 0.4,
                  "5,6": 0.3}
    result_add3 = is_additive(atpa_table)
    print("The additivity of ATPA synthase is", result_add3)


def is_ultrametric(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are ultrametric.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called. e.g. When comparing x and y, 
    also pass the threshold parameter: is_almost_equal(x, y, threshold).

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal
    Returns:
        (bool) True if the given distance metric is an ultrametric,
    False otherwise."""

    #
    # YOUR CODE HERE ()
    # construct new dictionary to store the changed keys and distances
    # dist is the original table, used_dict is the always changing table at every step,
    # dict stores all distances, including original ones and merges ones
    used_dict = dict = dist
    while len(used_dict) > 3:
        min_dist = 1
        # find the minimal distance
        for key in used_dict.keys():
            if used_dict[key] < min_dist:
                min_dist = used_dict[key]
        min_key = list(used_dict.keys())[list(used_dict.values()).index(min_dist)]
        # two sequences for the smallest distances
        first = min_key.split(',')[0]
        second = min_key.split(',')[-1]
        new_dict = {}
        # a new dictionary with merged sequences, first and second, deleted
        for key, value in used_dict.items():
            if first in key or second in key:
                continue
            else:
                new_dict[key] = value
        list_key = []
        # a list that stores all current sequences
        for i in used_dict.keys():
            key1 = i.split(',')[0]
            key2 = i.split(',')[1]
            # keys that contains sequences to be merged will not be included in the list
            if key1 not in list_key:
                list_key.append(key1)
            if key2 not in list_key:
                list_key.append(key2)
        # go through the sequence list and compute new distances
        for i in list_key:
            if str(i) in first or str(i) in second:
                continue
            else:
                # the new key is the combination of merged sequence and every other sequence in the list
                # the new key should let sequences stores according to the order
                if first < str(i):
                    new_key = first + '-' + second + ',' + str(i)
                elif first > str(i):
                    new_key = str(i) + ',' + first + '-' + second
                # when getting value from dict we should also be aware of the order of two sequence
                if str(i) > first:
                    score1 = dict[first + ',' + str(i)]
                elif first > str(i):
                    score1 = dict[str(i) + ',' + first]
                if str(i) > second:
                    score2 = dict[second + ',' + str(i)]
                elif second > str(i):
                    score2 = dict[str(i) + ',' + second]
                # the new value of new key is the average of distance between each sequence in merged sequence and the other sequence in list
                new_value = (score1 + score2)/2
                new_dict[new_key] = new_value
                dict[new_key] = new_value
        # update the current table
        used_dict = new_dict
    # check if one of the distances is smaller and the other two are equal
    if is_almost_equal(list(used_dict.values())[0],list(used_dict.values())[1],threshold) and list(used_dict.values())[0] > list(used_dict.values())[2]:
        return True
    elif is_almost_equal(list(used_dict.values())[1],list(used_dict.values())[2],threshold) and list(used_dict.values())[1] > list(used_dict.values())[0]:
        return True
    elif is_almost_equal(list(used_dict.values())[0],list(used_dict.values())[2],threshold) and list(used_dict.values())[0] > list(used_dict.values())[1]:
        return True
    return False


def is_additive(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are additive.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called. e.g. When comparing x and y, 
    also pass the threshold parameter: is_almost_equal(x, y, threshold).

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal

    Returns:
        (bool) Return True if the given distance metric is additive, 
        False otherwise."""

    #
    length = len(dist)
    a = 1
    b = -1
    c = - length * 2
    # this computes the number of sequences according to the length of input dictionary
    d = (-b + pow(b * b - 4 * a * c, 0.5)) / (2 * a)

    used_dict = dist
    avg_dist = {}  # averaged distances (r values) for all current nodes
    for i in range(1, int(d) + 1):
        avg_dist[str(i)] = 0
        for j in range(1, int(d) + 1):
            if i != j:
                avg_dist[str(i)] += dist["%d,%d" % (i, j)] if i < j else \
                    dist["%d,%d" % (j, i)]
        avg_dist[str(i)] = avg_dist[str(i)] / (len(list(used_dict.keys())) - 2)
    while d > 4:
        D = {}
        # Loop through each pair of nodes as entered in the dist dict
        # Use the entries from the avg_dist and calculate entries for the adjusted distance table D
        for i in used_dict.keys():
            first = i.split(',')[0]
            second = i.split(',')[1]
            new_value = used_dict[i] - avg_dist[first] - avg_dist[second]
            D[i] = new_value
        minimal = 0
        i = ''
        j = ''
        # find the minimal adjusted distance
        for key, value in D.items():
            if value < minimal:
                minimal = value
                # i and j are two sequences of the key
                i = key.split(',')[0]
                j = key.split(',')[1]
        # construct new table with sequence that will be merged excluded
        new_dict = {}
        for key, value in used_dict.items():
            if i in key or j in key:
                continue
            else:
                new_dict[key] = value
        # build a list that stores all current sequences
        list_key = []
        for n in used_dict.keys():
            key1 = n.split(',')[0]
            key2 = n.split(',')[1]
            if key1 not in list_key:
                list_key.append(key1)
            if key2 not in list_key:
                list_key.append(key2)
        # k is the merged sequence
        k = i + '-' + j
        # go through the sequence list and compute new distances
        for m in list_key:
            if m != i and m != j:
                # the new key is the combination of merged sequence and every other sequence in the list
                # the new key should let sequences stores according to the order
                dist_ij = dist[str(i) + ',' + str(j)]
                if m < k:
                    new_key = m + ',' + str(k)
                elif m > k:
                    new_key = str(k) + ',' + m
                # when getting values from dict we should also be aware of the order of two sequences
                if i > m:
                    dist_im = dist[m + ',' + str(i)]
                elif m > i:
                    dist_im = dist[str(i) + ',' + m]
                if j > m:
                    dist_jm = dist[m + ',' + str(j)]
                elif m > j:
                    dist_jm = dist[str(j) + ',' + m]
                new_value = (dist_im + dist_jm - dist_ij) / 2
                new_dict[new_key] = new_value
                dist[new_key] = new_value
        used_dict = new_dict
        temp = list_key
        temp.remove(i)
        temp.remove(j)
        # node_list is the list for all sequences,
        # excluding the two original sequences for merging, including the new merged sequence
        node_list = temp
        node_list.append(k)
        # update the average dictionary
        if len(node_list) > 2:
            avg_dist[k] = 0
            for ind in range(0, len(node_list) - 1):
                m = node_list[ind]
                avg_dist[m] = avg_dist[m] * (len(node_list) - 1)
                avg_dist[m] -= dist["%s,%s" % (m, i)] if m < i \
                    else dist["%s,%s" % (i, m)]
                avg_dist[m] -= dist["%s,%s" % (m, j)] if m < j \
                    else dist["%s,%s" % (j, m)]
                avg_dist[m] += dist["%s,%s" % (m, k)] if m < k \
                    else dist["%s,%s" % (k, m)]

                avg_dist[m] /= (len(node_list) - 2)
                avg_dist[k] += dist["%s,%s" % (m, k)] if m < k \
                    else dist["%s,%s" % (k, m)]

            avg_dist[k] = avg_dist[k] / (len(node_list) - 2)
        d = len(node_list)
    sums = []
    # check the sums of scores of all combination of sequences
    for i in used_dict.keys():
        key1 = i.split(',')[0]
        key2 = i.split(',')[1]
        for j in used_dict.keys():
            if i != j and key1 not in j and key2 not in j and i < j:  # avoid repetition
                sums.append(used_dict[i] + used_dict[j])
    # check if two of the sums are equal and the other one is smaller
    if is_almost_equal(sums[0],sums[1],threshold) and sums[0] > sums[2]:
        return True
    elif is_almost_equal(sums[1],sums[2],threshold) and sums[1] > sums[0]:
        return True
    elif is_almost_equal(sums[0],sums[2],threshold) and sums[0] > sums[1]:
        return True
    return False


def is_almost_equal(num_1, num_2, threshold):
    """
    Return true if the difference between the two parameters is negligible
    enough that the parameters can be considered equal.
    """
    return abs(num_1 - num_2) <= threshold


if __name__ == '__main__':
    solve_ultrametric_additive()
