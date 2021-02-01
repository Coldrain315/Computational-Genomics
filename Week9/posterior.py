import sys
from math import log, exp
from compsci260lib import get_fasta_dict


def run_posterior():
    # input_file = "artificial.genome.fasta"
    input_file = "bacterial.genome.fasta"
    hmm_file = "HMM.methanococcus.txt"
    posterior = posterior_decoding(input_file, hmm_file)

    # Report the first and last ten segments in your decoding
    #
    print("The first ten segments:")
    for i in range(10):
        print("The start position is %d, the end position is %d, and the state is %s" % (
            posterior[i]['start'], posterior[i]['end'], posterior[i]['state']))
        # print("The last ten segments are ", posterior[-10:])
    print("The last ten segments:")
    for i in range(9, -1, -1):
        j = len(posterior) - i - 1
        print("The start position is %d, the end position is %d, and the state is %s" % (
            posterior[j]['start'], posterior[j]['end'], posterior[j]['state']))
        # print("The last ten segments are ", posterior[-10:])
    counts = count_segments(posterior)
    print("The number of segments that exist in state1 is %d, the number of segments that exist in state2 is %d" % (
        counts['state1'], counts['state2']))

    #


def posterior_decoding(input_file, hmm_file):
    """
    Calculate the posterior decoding and return the decoded segments.

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

    # Read in the input files
    f_hmm_file = open(hmm_file)
    if f_hmm_file is None: sys.exit("Can't open file: " + hmm_file)

    # read the state names and number
    states = f_hmm_file.readline().split()
    K = len(states)

    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]

    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row

    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()

    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row

    f_hmm_file.close()

    seq_dict = get_fasta_dict(input_file)
    emit_str = list(seq_dict.values())[0]  # there's only 1

    print("Done reading sequence of length ", len(emit_str))

    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emitted_symbols,
                          emit_probs, emit_str)

    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions,
                            emitted_symbols, emit_probs, emit_str)

    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        for k in range(K):
            posterior[i][k] = forward[i][k] + backward[i][k]
    # Create the list of decoded segments to return
    #
    list_dict = []
    start = 1
    max_prob = [posterior[0][k] for k in range(K)]
    prev_state = max_prob.index(max(max_prob))
    for i in range(1, len(emit_str)):
        max_prob = [posterior[i][k] for k in range(K)]
        state = max_prob.index(max(max_prob))
        if state != prev_state:
            end = i
            list_dict.append({'start': start, 'end': end, 'state': 'state' + str(state + 1)})
            start = i + 1
            prev_state = state
    #
    return list_dict


def run_forward(states, initial_probs, transitions, emitted_symbols,
                emit_probs, emit_str):
    """Calculates the forward probability matrix.

    Arguments:
        states (list of str): list of states as strings
        initial_probs (list of float): list of initial probabilities for each
            state
        transitions (list of list of float): matrix of transition probabilities
        emission_symbols (list of str): list of emission symbols
        emit_probs (list of list of float): matrix of emission probabilities
            for each state and emission symbol
        emit_str (str):

    Returns:
        (list of list of floats): matrix of forward probabilities
    """

    K = len(states)
    N = len(emit_str)
    forward = [[float(0) for _ in range(K)] for _ in range(N)]

    # Initialize
    emit_index = get_emit_index(emit_str[0], emitted_symbols)
    for k in range(K):
        forward[0][k] = log(initial_probs[k]) + log(emit_probs[k][emit_index])
    # print(forward)
    # Iterate
    for i in range(1, N):
        emit_index = get_emit_index(emit_str[i].upper(), emitted_symbols)
        for a in range(K):
            for b in range(K):  # Compute the forward probabilities for the states
                prob = forward[i - 1][b] + log(transitions[b][a]) + log(emit_probs[a][emit_index])
                if forward[i][a] == 0:
                    forward[i][a] = log_sum([prob])
                else:
                    forward[i][a] = log_sum([forward[i][a], prob])

    return forward


def run_backward(states, initial_probs, transitions, emitted_symbols,
                 emit_probs, emit_str):
    """Calculates the backward probability matrix.

        Arguments:
            states (list of str): list of states as strings
            initial_probs (list of float): list of initial probabilities for 
                each state
            transitions (list of list of float): matrix of transition
                probabilities
            emission_symbols (list of str): list of emission symbols
            emit_probs (list of list of float): matrix of emission
                probabilities for each state and emission symbol
            emit_str (str):

        Returns:
            (list of list of floats): matrix of backwards probabilities
    """

    K = len(states)
    N = len(emit_str)

    backward = [[float(0) for _ in range(K)] for _ in range(N)]

    # Initialize
    for k in range(K):
        backward[N - 1][k] = log(1)  # which is zero, but just to be explicit...

    # Iterate
    # Compute the backward probabilities for the states
    for i in range(N - 2, -1, -1):
        emit_index = get_emit_index(emit_str[i + 1].upper(), emitted_symbols)
        for a in range(K):
            for b in range(K):
                prob = backward[i + 1][b] + log(transitions[a][b]) + log(emit_probs[b][emit_index])
                if backward[i][a] == 0:
                    backward[i][a] = log_sum([prob])
                else:
                    backward[i][a] = log_sum([backward[i][a], prob])
    return backward


def log_sum(factors):
    """
    Calculates the sum for log values
    Arguments: list of factors to be sum up
    Returns: the log sum of factors
    """
    if abs(min(factors)) > abs(max(factors)):
        a = min(factors)
    else:
        a = max(factors)
    total = 0
    for x in factors:
        total += exp(x - a)
    return a + log(total)


def get_emit_index(input_val, alphabet):
    """does a linear search to find the index of a character"""
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    sys.exit("Could not find character " + input_val)


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


if __name__ == '__main__':
    """Call run_posterior(), do not modify"""
    run_posterior()
