
from compsci260lib import *
import sys, random, os


def generate_HMM_sequence(hmm_file, seq_length):
    """
    Load an HMM specification from file, and generate state and observed
    sequences through sampling of the HMM.

    The HMM will have states labeled as strings and observations as characters
    which will be useful when we generate an HMM with nucleotide sequence
    observations.

    An example return sequence may look like:
        (['W', 'S', 'W'], ['A', 'T', 'G']) or
        (['F', 'L', 'L'], ['2', '6', '6'])

    Arguments:
        hmm_file (str): path to the HMM file
        seq_length (int): the length of the sequence we will generate

    Returns:
        a tuple of
        (the state sequence as strings,
         observed sequence as single character strings)
    """

    if not os.path.exists(hmm_file):
        print("Can't open HMM parameter file: %s" % hmm_file)
        return -1

    f = open(hmm_file, "r")

    # read the state names
    states = f.readline().strip().split()

    # read the initial probabilities
    initial_probs = f.readline().strip().split()

    initial_probs = [float(p) for p in initial_probs]
    # read the transition matrix
    transitions = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        transitions[state] = [float(p) for p in matrix_row]
    # read the input alphabet
    input_alphabet = f.readline().strip().split()
    # read the emission matrix
    emission = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        emission[state] = [float(p) for p in matrix_row]
        # normalize
        sum_emission = sum(emission[state])
        emission[state] = [x / sum_emission for x in emission[state]]
    f.close()
    #
    state_seq = []
    obs_seq = []
    new_trans_rate = []
    new_emis_rate = {}
    # the range of prob is 0 < prob < 1, however, the range of rate used for function random index is 0-100
    # thus we conduct the multiplication first
    for i, prob in enumerate(initial_probs):
        new_trans_rate.append(prob * 100)
    for state, lists in emission.items():
        new_emis_rate[state] = []
        for prob in lists:
            new_emis_rate[state].append(prob * 100)
    # initialize the first state and first ob
    initial_state = states[random_index(new_trans_rate)]
    initial_obs = input_alphabet[random_index(new_emis_rate[initial_state])]
    # add the first state and first ob to sequences
    state_seq.append(initial_state)
    obs_seq.append(initial_obs)

    cur_state = initial_state
    # go through the process and generate sequences according to the required sequence length
    for i in range(seq_length-1):
        # change the rate range to 100
        cur_trans_rate = [i*100 for i in transitions[cur_state]]
        # randomly get current state based on the probability
        cur_state = states[random_index(cur_trans_rate)]
        state_seq.append(cur_state)
        # randomly get current observation based on the probability
        cur_obs = input_alphabet[random_index(new_emis_rate[cur_state])]
        obs_seq.append(cur_obs)
    return (state_seq, obs_seq)  # a tuple containing the state and observation sequences


def random_index(rate):
    """generate random variable according to the given probability
       Argument: a list of probability of events, in range 100; the order are accordant to the list of events
       Return: the index of the randomly chosen event"""
    start = 0
    index = 0
    randnum = random.randint(1, sum(rate))
    for index, scope in enumerate(rate):
        start += scope
        if randnum <= start:
            break
    return index


print(generate_HMM_sequence("HMM.parameters.txt",4))
