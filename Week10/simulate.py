import sys
import random
import numpy as np
from compsci260lib import *

def run_simulate():
    """
    Simulates the sequencing process then empirically compute and report 
    the quantities from parts a-c.
    """

    iterations = 20
    G = 3*1000000
    R = 4*10000
    L = 450

    # Call simulate(G, R, L) `iteration` times. Report the empirical values 
    # from parts a-c for each iteration
    #
    coverage = 0
    unsequenced = 0
    n_contig = 0
    len_contig = 0
    for i in range(iterations):
        result = simulate(G, R, L)
        coverage += result[0]
        unsequenced += result[1]
        n_contig += result[2]
        len_contig += result[3]
        print("The coverage C is %f, the number of nucleotides in the genome that remain unsequenced is %d,"
              "the number of contigs is %d, the length of each contig is %f." % result)
    print("The average coverage C is %f, the average number of nucleotides in the genome that remain unsequenced is %f,"
          "the average number of contigs is %f, the average length of each contig is %f" % (coverage/20, unsequenced/20, n_contig/20, len_contig/20))
    #


def simulate(G, R, L):
    """
    Simulates one iteration of the sequencing process and empirically compute 
    the empirical coverage (average number of times a nucleotide in the genome
    was sequenced), the number of nucleotides not covered by any read, the 
    number of contigs assuming you can use the oracular assembly algorithm to
    assemble all the reads, and the average length of these contigs.

    Args:
        G (int) - the length of the genome
        R (int) - the number of reads
        L (int) - the length of each read

    Returns
        a tuple of floats:

            (Empirical coverage, 
             Number of nucleotides not covered by any read, 
             Number of contigs,
             Average length of these contigs)
    """

    # initialize the array for storing number of times of covering
    time_list = np.zeros(G+L)

    # randomly pick R reads with length L
    for _ in range(R):
        selected_index = random.randint(0,G+L)
        time_list[selected_index:selected_index+L] += 1
    # deal with edge effect
    time_list = time_list[L:G+L]
    average = sum(time_list)/G  # average number of each nucleotide's being covered, also known as coverage
    none = 0
    length = 0  # length of each contig
    num_contig = []  # list that stores the length of each contig, the length of num_contig is the number of contigs

    for i in time_list:
        if i == 0:  # 0 means not covered
            none += 1  # count the number of nucleotides not covered
            if length != 0:  # this means it is the turning point from contig to not covered part
                num_contig.append(length)  # record the length of the contig that ends here
            length = 0  # if not covered, there is no contig at current location, thus 0 length for contig
        else:
            length += 1  # not 0 means current location is covered, thus increment the length of current contig
    num_contig.append(length)
    #
    return (average, none, len(num_contig), sum(num_contig)/len(num_contig))


if __name__ == '__main__':
    """Call run_simulate(), do not modify"""
    run_simulate()
