from compsci260lib import *
import timeit
import random

# random_list needs to be available as a global variable
# for the timeit function to work properly
random_list = None


def brute_force(arr):
    """Get the maximum similarity score using brute force.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """
    # sum_list = []
    maximum = 0
    for i in range(len(arr)):
        temp_sum = 0
        for j in range(i, len(arr)):
            temp_sum += arr[j]
            if temp_sum > maximum:
                maximum = temp_sum
    return maximum


def divide_conquer(arr):
    """Get the maximum similarity score using divide and conquer.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """

    """Your code goes here..."""
    maximum = maxSubArraySum(arr, 0, len(arr) - 1)
    return maximum


# this function aims to find the max crossing sub-array
def maxCrossingSum(arr, l, m, h):
    # Include elements on left of mid.
    sm = 0
    left_sum = -float('inf')
    for i in range(m, l - 1, -1):
        sm = sm + arr[i]
        if sm > left_sum:
            left_sum = sm
            # Include elements on right of mid
    sm = 0
    right_sum = -float('inf')
    for i in range(m + 1, h + 1):
        sm = sm + arr[i]
        if sm > right_sum:
            right_sum = sm
            # Return sum of elements on left and right of mid
    # returning only left_sum + right_sum will fail for [-2, 1]
    return max(left_sum + right_sum, left_sum, right_sum)


# Returns sum of maximum sum sub-array in aa[l..h]
def maxSubArraySum(arr, l, h):
    # Base Case: Only one element
    if l == h:
        return arr[l]
    # Find middle point
    m = (l + h) // 2
    # Return maximum of following three possible cases
    return max(maxSubArraySum(arr, l, m),  # Maximum sub-array sum in left half
               maxSubArraySum(arr, m + 1, h),  # Maximum sub-array sum in right half
               maxCrossingSum(arr, l, m, h))  # c) Maximum sub-array sum such that the sub-array crosses the midpoint


def linear(arr):
    """Get the maximum similarity score in linear time.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """

    """Your code goes here..."""
    max_including_here = 0
    max_so_far = 0
    for a in arr:
        max_including_here = max(0, max_including_here + a)
        max_so_far = max(max_so_far, max_including_here)
    return max_so_far


def run_collective_similarity():
    """
    Run collective similarity and collect timing for each algorithm.

    Note: there is no reporting in this problem, `brute_force`, 
    `divide_conquer`, and `linear` will be evaluated for correctness.
    """

    # You can use this to test the correctness of your code by using
    # sample_list as an input to each function.

    # declare random_list as the same global variable defined
    # on line 7 for use in the timeit library
    global random_list

    sample_list = [2, -3, -4, 4, 8, -2, -1, 1, 10, -5]

    print("The maximum sum of sub-list found by brute force algorithm is", brute_force(sample_list))
    print("The maximum sum of sub-list found by divide and conquer algorithm is", divide_conquer(sample_list))
    print("The maximum sum of sub-list found by linear algorithm is", linear(sample_list))

    # This part below is used to test the runtime of your code, an example is
    # given below for brute force algorithm with a random list of length 100.
    # You will have to measure the runtime of each algorithm on every input size
    # given in the problem set.

    # """
    allowed_scores = [i for i in range(-10, 11)]

    for i in range(2, 6):
        length = 10 ** i
        random_list = [random.choice(allowed_scores) for x in range(length)]
        bruteforce_runtime = timeit.timeit('brute_force(random_list)',
                                           setup="from __main__ import brute_force, random_list", number=1)
        print("For the brute force algorithm, its running time on sets of random inputs of length 10^%d is %4f" % (
            i, bruteforce_runtime))

    for i in range(2, 8):
        length = 10 ** i
        random_list = [random.choice(allowed_scores) for x in range(length)]
        dq_runtime = timeit.timeit('divide_conquer(random_list)',
                                   setup="from __main__ import divide_conquer, maxCrossingSum, maxSubArraySum, random_list",
                                   number=1)
        print("For the divide and conquer algorithm, its running time on sets of random inputs of length 10^%d is %4f" % (
                i, dq_runtime))

    for i in range(2, 9):
        length = 10 ** i
        random_list = [random.choice(allowed_scores) for x in range(length)]
        linear_runtime = timeit.timeit('linear(random_list)',
                                       setup="from __main__ import linear, random_list",
                                       number=1)
        print("For the linear algorithm, its running time on sets of random inputs of length 10^%d is %4f" % (
            i, linear_runtime))


if __name__ == '__main__':
    """Run run_collective_similarity(). Do not modify this code"""
    run_collective_similarity()
