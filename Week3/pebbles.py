from compsci260lib import *


def solve_pebble(grid_file):
    """Code for the "Pebble Beach" problem. This problem involves implementing
    an O(n) dynamic programming algorithm for computing the maximum value of
    the placement of pebbles under the constraint that no pebbles can be
    vertically or horizontally adjacent.

    Args: grid_file (str): a string with the name of the file that contains
          the grid of values. Each line of that file should contain a row
          of four integers, separated by tabs

    Returns: the maximal score for the optimal pebble placements
    """

    """Your code goes here..."""
    grid_list = read_file(grid_file)
    # find all possible combinations; 1 means a pebble here, 0 is null
    all_patterns = [[a, b, c, d] for a in [0, 1] for b in [0, 1] for c in [0, 1] for d in [0, 1]]
    possible_patterns = []
    count_conditions = 0
    for a, b, c, d in all_patterns:
        if a + b != 2 and b + c != 2 and c + d != 2:  # not adjacent
            possible_patterns.append([a, b, c, d])
            count_conditions += 1

    pattern_type_match = {}
    for i in range(count_conditions):
        temp = []
        for j in range(count_conditions):
            if examine(possible_patterns[i], possible_patterns[j]):
                temp.append(j)
        pattern_type_match[i] = temp  # pattern_type_match is a dictionary that
        # contains all possible pattern that match as the next row for each pattern type

    possible_sums = [[0 for i in range(count_conditions)]]
    for i in range(len(grid_list)):
        temp_sum = []
        for j in possible_patterns:
            sums = 0
            for index, p in enumerate(j):
                sums += int(p) * int(grid_list[i][index])
            temp_sum.append(sums)
        possible_sums.append(temp_sum)
    dp = possible_sums  # dp will be the 2d-array that stores current max

    for i in range(1, len(grid_list) + 1):
        ans = 0
        for j in pattern_type_match.keys():  # stands for the pattern type of current row
            mx = 0
            for k in pattern_type_match[j]:  # stands for the pattern type of previous row that can match the current row
                mx = max(mx, possible_sums[i - 1][k])  # takes the maximum of previous row
                # the max of current row is determined by its previous row
            dp[i][j] = mx + possible_sums[i][j]  # sum of current position and max previous row
            ans = max(ans, dp[i][j])  # answer is the maximum sum of all positions that has been gone through
    return ans


def read_file(grid_file):
    lines = open(grid_file).read().split("\n")
    grids = []
    for line in lines:  # reads lists of data from txt line by line
        grid = line.split("\t")
        if len(grid) < 4:
            continue
        else:
            grids.append(grid)
    return grids


def examine(condition1, condition2):
    for index, send in enumerate(condition1):
        if send*condition2[index] == 1:  # this means the
            return False
    return True


def run_solve_pebble():
    """
    Run solve pebbles, you may try an create different grid files to debug
    that match the formatting of grid.txt (tab separated values)

    Note: there is no reporting in this problem, only `solve_pebbles` will
    be evaluated for correctness.
    """
    max_score = solve_pebble('grid.txt')
    print("The max score for this grid is %d" % max_score)


if __name__ == '__main__':
    """Run run_solve_pebble(). Do not modify this code"""
    run_solve_pebble()
