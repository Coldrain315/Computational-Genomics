from compsci260lib import *
from generate_HMM_sequence import generate_HMM_sequence


def run_analyze_sequence():
    """Load the nucleotide sequence HMM file and generate a 100000 length sequence,
    then do some analysis of the results, including the average length in each state"""
    seq_length = 100000
    (state_sequence,
     observed_sequence) = generate_HMM_sequence("HMM.parameters.txt",
                                                seq_length)
    state_str = ''
    observed_str = ''
    for i in state_sequence:
        state_str += i
    for i in observed_sequence:
        observed_str += i

    # Write your generated sequence to file.
    filename = "nucleotide_sequence.txt"
    f = open(filename, "w")

    if f is None:
        sys.stderr("output file not found: " + filename)
    else:
        # output the results
        f.write(state_str + " \n")
        f.write(observed_str + " \n")
        f.close()
        print("Wrote %s" % filename)

    # Report the states that prevailed whenever
    # TCGA was observed, and their frequencies when observing TCGA.
    subs = ['T', 'C', 'G', 'A']
    frequencies = compute_state_frequencies(state_sequence, observed_sequence, subs)
    print(frequencies)
    for key, value in frequencies.items():
        print("The frequency for sequences %s when observing TCGA is %d" % (key, value))
    # Finally, report the average length spent in each state across the entire sequence.
    #
    length = compute_average_length(state_sequence)
    for key, value in length.items():
        print("The average length spent in %s is %d. " % (key, value))

    #


def compute_average_length(state_sequence):
    """
    Given a state sequence, return the average length in each state

    Arguments:
        state_sequence (list of one char strings): generated sequence of states

    Returns:
        a dictionary mapping the state name to the average length in the state.
        Example return:
        {
            'W': 5.5,
            'S': 4.5
        }

    """

    dictionary = {}  # the returning dictionary mapping the state name to the average length in the state.
    temp_dict = {}  # the temporary dictionary storing frequency and segment number of state name
    i = 1

    while True:
        while i < len(state_sequence):
            cur_char = state_sequence[i - 1]
            next_char = state_sequence[i]
            if cur_char == next_char:  # if current char is same as next char, continue to go to next char; loop until a difference is met, and then break
                i += 1
                if cur_char not in temp_dict.keys():
                    temp_dict[cur_char] = {'long': 1, 'count': 0}
                else:
                    temp_dict[cur_char]['long'] += 1
            else:  # if current char is different from next char, the state is changed
                if cur_char not in temp_dict.keys():  # if the first time the key appears, initialize the key and its values
                    temp_dict[cur_char] = {'long': 1, 'count': 0}
                else:  # if the current substring already exist as key, add 1 to its count
                    temp_dict[cur_char]['long'] += 1
                i += 1
                break
        temp_dict[cur_char]['count'] += 1  # every time a loop is break, it means the state change, thus add one to the count of segment of the previous state
        if i >= len(state_sequence):
            break
    # there will be one final char left because the loop can not go through the end one
    # similar process of key initialization and adding
    if state_sequence[-1] == state_sequence[-2]:
        temp_dict[cur_char]['long'] += 1
    elif state_sequence[-1] != state_sequence[-2]:
        if state_sequence[-1] not in temp_dict.keys():
            temp_dict[state_sequence[-1]] = {'long': 1, 'count': 0}
        else:
            temp_dict[state_sequence[-1]]['long'] += 1
            temp_dict[state_sequence[-1]]['count'] += 1
    # go through the temporary and calculate average length by dividing total length by number of segments appeared
    for key, value in temp_dict.items():
        avrg = value['long'] / value['count']
        dictionary[key] = avrg
    return dictionary
    #


def compute_state_frequencies(st_sequence, observed_sequence, subsequence):
    """Given the state and observed sequences, return the state sequences that
    emitted the query subsequence and frequency in which those state sequences
    sequences emitted the subsequence.

    Arguments:
        state_sequence (list of one char strings): generated sequence of states
        observed_sequencs (list of one char strings): generated observed 
            sequence
        subsequence (list of one char strings): the observed subsequence to
            count

    Returns:
        a dictionary mapping the state name to the frequency of observing the
        provided sequence. Example return:
        {
            'WWWW': 2,
            'WWWS': 1,
            ...
        }
    """

    #
    index = 0
    dictionary = {}
    # loop through all observed sequences
    while index < len(observed_sequence):
        count = 0
        # loop through the subsequence
        for sub in subsequence:
            if sub == observed_sequence[index]:  # compare if current sequence in observed sequence is same as the current one in subsequence
                count += 1  # add one to the count of matched sequence
                index += 1  # continue to next sequence of observed sequence
                if index >= len(observed_sequence):
                    break
            else:  # break if current sequence in observed sequence does not match the current one in subsequence
                break
        if count == len(subsequence):  # if all sequences in the given subsequence have been matched
            state = ''
            for i in st_sequence[index - len(subsequence):index]:
                state += i
            if state not in dictionary.keys():
                dictionary[state] = 1
            else:
                dictionary[state] += 1
        else:  # if not matched, go back to the point before comparing to subsequence one by one and continue to move to next char in observed sequence
            index = index-count+1
    return dictionary
    #


if __name__ == '__main__':
    """Main method call, do not modify"""
    run_analyze_sequence()
