import sys
from math import log
from compsci260lib import *
import math


def run_viterbi():
    hmm_file = 'HMM.methanococcus.txt'
    input_file = 'bacterial.genome.fasta'
    print("Decoding the Viterbi path for %s using %s" % (input_file, hmm_file))

    vit_ret = viterbi_decoding(input_file, hmm_file)

    # vit_ret = viterbi_decoding("try.genome.fasta", "HMM.try.txt")
    # Collect the segment counts for each state
    counts = count_segments(vit_ret)
    # Report the first 10 and last 10 segments of your decoding
    #
    # print("The first ten segments are ", vit_ret[0:10])
    # print("The last ten segments are ", vit_ret[-10:])
    print("The first ten segments:")
    for i in range(10):
        print("The start position is %d, the end position is %d, and the state is %s" % (
            vit_ret[i]['start'], vit_ret[i]['end'], vit_ret[i]['state']))
        # print("The last ten segments are ", posterior[-10:])
    print("The last ten segments:")
    for i in range(9, -1, -1):
        j = len(vit_ret) - i - 1
        print("The start position is %d, the end position is %d, and the state is %s" % (
            vit_ret[j]['start'], vit_ret[j]['end'], vit_ret[j]['state']))
    #

    # Then report the number of segments that exist in each state.
    #
    print("The number of segments that exist in state1 is %d, the number of segments that exist in state2 is %d" % (counts['state1'], counts['state2']))
    #


def viterbi_decoding(input_file, hmm_file):
    """Calculate the viterbi decoding of an input sequence

    Arguments:
        input_file (str): path to input fasta file
        hmm_file (str): path to HMM file

    Returns:
        A list of dictionaries of segments in each state (1-indexed).
        An example output may look like:

        [
            {‘start’: 1, ‘end’: 12, ‘state’: ‘state2’},
            {‘start’: 13, ‘end’: 20, ‘state’: ‘state1’},
            ...
        ]
    """

    # open hmm file
    try:
        f_hmm_file = open(hmm_file, 'r')
    except IOError:
        print("IOError: Unable to open HMM file: %s." % hmm_file)
        print("Exiting.")
        sys.exit(1)

    # read the state names
    states = f_hmm_file.readline().split()
    K = len(states)

    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]  # here is 0.5 and 0.5 in methanococcus

    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row
    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()
    # print(emitted_symbols)
    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row
    # print("emit_prob", emit_probs)
    f_hmm_file.close()

    seq_dict = get_fasta_dict(input_file)
    # the observation sequence
    emit_str = list(seq_dict.values())[0]  # there's only 1

    print("Read a sequence of length", len(emit_str))

    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]  # K rows and len(emit_str) columns
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)]

    # Initialize the first column of the matrix
    for i in range(K):
        in_index = get_emit_index(emit_str[0].upper(), emitted_symbols)
        viterbi[i][0] = log(emit_probs[i][in_index]) + log(initial_probs[i])
        pointers[i][0] = str(i + 1) + '-' + str(i + 1)
        # print(i, in_index, emit_probs[i][in_index], initial_probs[i])
        # print(viterbi[i][0])
    # Build the matrix column by column
    for j in range(1, len(emit_str)):
        in_index = get_emit_index(emit_str[j].upper(), emitted_symbols)
        prev_viterbi = [viterbi[n][j - 1] for n in range(K)]
        for i in range(K):  # the current state
            # Compute the entries viterbi[i][j] and pointers[i][j]
            # Tip: Use float('-inf') for the value of negative infinity
            #
            prev_trans = [log(transitions[n][i]) for n in range(K)]  #
            prev_sum = [prev_trans[n] + prev_viterbi[n] for n in range(K)]
            max_prev = max(prev_sum)
            prev_state = prev_sum.index(max_prev)
            viterbi[i][j] = log(emit_probs[i][in_index]) + max_prev
            pointers[i][j] = str(prev_state + 1) + '-' + str(i + 1)
            # 
            pass  # placeholder line
    # print(viterbi, pointers)
    max_prob = [viterbi[i][len(emit_str) - 1] for i in range(K)]
    state = max_prob.index(max(max_prob))
    # Traceback, stored as a list segment lengths in each state as dictionaries
    nth = len(emit_str)
    seg_dicts = []
    end = nth
    while nth >= 1:
        pointer = pointers[state][nth - 1]
        prev_state = pointer.split('-')[0]
        cur_state = pointer.split('-')[1]
        if prev_state != cur_state or nth == 1:
            seg_dicts.append({'start': nth, 'end': end, 'state': "state" + cur_state})
            end = nth - 1
        state = int(prev_state) - 1
        nth -= 1
    #
    # print(sorted(seg_dicts, key=lambda keys: keys['start']))
    return sorted(seg_dicts, key=lambda keys: keys['start'])


def count_segments(vit_ret):
    """Calculate the number of segments appearing in each state of
    the viterbi path

    Arguments:
        vit_ret (list of dicts): dictionary of lengths in each state.
            see: return value of viterbi_decoding

    Returns:
        a dictionary of states to number of occurrences in the state. 
        e.g. {'state1': 10, 'state2': 9}
    """

    #
    occur_dict = {'state1': 0, 'state2': 0}
    length_dict = {'state1': 0, 'state2': 0}
    for i, dic in enumerate(vit_ret):
        if int(dic['state'][-1]) == 1:
            occur_dict['state1'] += 1
            length_dict['state1'] += (dic['end'] - dic['start'] + 1)
        else:
            occur_dict['state2'] += 1
            length_dict['state2'] += (dic['end'] - dic['start'] + 1)
    #
    return occur_dict


def get_emit_index(input_val, alphabet):
    """Get the index of the emission value in the alphabet

    Note: This will be a useful function for indexing the emission
    probabilities"""
    return alphabet.index(input_val)


if __name__ == '__main__':
    """Call run_viterbi(), do not modify"""
    run_viterbi()
