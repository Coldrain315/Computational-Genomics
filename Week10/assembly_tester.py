from compsci260lib import *


def run_assembly_tester():
    """
    Read in the contig files, verify the reads in each contig, and estimate
    the gap length between the contigs.
    """
    reads_file = 'paired.reads.fasta'
    contig0_file = 'contig0.fasta'
    contig1_file = 'contig1.fasta'

    # load the fasta files
    reads = get_fasta_dict(reads_file)
    contig0 = get_fasta_dict(contig0_file)
    contig1 = get_fasta_dict(contig1_file)

    # Determine how the reads map to the contigs
    contig_reads = find_contig_reads(reads, contig0, contig1)

    # Report the reads whose ends map to both contig0 and contig1
    # and their estimated gap lengths
    #
    # d, e)
    count = 0
    for i in contig_reads:
        temp_dicta = contig_reads[i]['contig_a']
        temp_dictb = contig_reads[i]['contig_b']
        if temp_dicta != temp_dictb:
            count += 1
        # elif temp_dicta is None and temp_dictb is None:
        #     print("The sequence %s doesn’t belong to the assembly." % i)
    print("There are %d pair of reads that map to both contig0 and contig1." % count)

    # supercontig
    supercontig = get_fasta_dict("supercontig.fasta")['super_contig0']
    con0 = supercontig.split("<---bridged gap--->")[0]
    con1 = supercontig.split("<---bridged gap--->")[1]
    seq = con0 + con1
    for key, value in contig_reads.items():
        if value['contig_a'] != value['contig_b']:
            start = seq.index(reads[key+"a"]) + 1
            end = seq.index(reads[key+"b"]) + len(reads[key+"b"])
            length = end - start + 1
            print("The length of the sequence %s is %d, thus the gap range is %d to %d" %(key,length,1990-length,2010-length))


def find_contig_reads(reads, contig0, contig1):
    """
    Determine whether the sequencing reads map to contig0/contig1/both/neither
    and where in the contig it matches. `reads` will contain both ends 'a' and
    'b', but you will return a dictionary using the read name as a whole
    (without 'a' and 'b').

    It should return a dictionary mapping the name of the mated pair of reads
    to:
        - 'contig_a' (str): the contig in which read 'a' was found, as 
          `contig0', `contig1', or None.
        - start_a (int): the start position (1-indexed) read end 'a' mapped to
          within its respective contig (None if not found in any contig)
        - end_a (int): the end position (1-indexed) read end 'a' mapped to
          within its respective contig (None if not found in any contig)
        - 'contig_b' (str): the contig in which read 'b' was found, as
          `contig0', `contig1', or None.
        - start_b (int): the start position (1-indexed) read end 'b' mapped to
          within its respective contig (None if not found in any contig)
        - end_b (int): the end position (1-indexed) read end 'b' mapped to
          within its respective contig (None if not found in any contig)

    The returned dictionary should look something like:
    {
        'seq1': {
            'contig_a': 'contig0',
            'start_a': 301,
            'end_a': 800,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        'seq2': {
            'contig_a': 'contig1',
            'start_a': 1101,
            'end_a': 1600,
            'contig_b': 'contig0'
            'start_b': 201,
            'end_b': 700,
        },
        'seq3' : {
            'contig_a': None,
            'start_a': None,
            'end_a': None,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        ...
    }

    Arguments:
        reads (dict str to str): dictionary of sequencing reads
        contig0 (dict str to str): dictionary of reads in contig0
        contig1 (dict str to str): dictionary of reads in contig1

        see: get_fasta_dict

    Returns:
        Dictionary mapping reads to information about their locations in contigs.
    """

    #
    dictionary = {}
    for i in reads:
        num = i[0:-1]
        if num not in dictionary.keys():
            dictionary[num] = {}
        sequence = reads[i]
        ab = "contig_" + i[-1]
        start = "start_" + i[-1]
        end = "end_" + i[-1]
        if sequence in contig0['contig0']:
            dictionary[num][ab] = 'contig0'
            dictionary[num][start] = contig0['contig0'].index(sequence) + 1
            dictionary[num][end] = contig0['contig0'].index(sequence) + len(sequence)
        elif sequence in contig1['contig1']:
            dictionary[num][ab] = 'contig1'
            dictionary[num][start] = contig1['contig1'].index(sequence) + 1
            dictionary[num][end] = contig1['contig1'].index(sequence) + len(sequence)
        else:
            dictionary[num][ab] = None
            dictionary[num][start] = None
            dictionary[num][end] = None
    for i in dictionary:
        temp_dicta = dictionary[i]['contig_a']
        temp_dictb = dictionary[i]['contig_b']
        if temp_dicta is None and temp_dictb is None:  # doesn't belong  to the assembly
            print("The sequence %s doesn’t belong to the assembly." % i)
        elif temp_dictb == temp_dicta:  # appear in the same contig
            if dictionary[i]['end_b'] - dictionary[i]['start_a'] + 1 < 1990 or dictionary[i]['end_b'] - dictionary[i]['start_a'] + 1> 2010:
                print("The sequence %s doesn’t obey the second rule." % i)
        elif temp_dictb != temp_dicta:
            if temp_dicta == contig1 and temp_dictb == contig0:
                print("The sequence %s doesn’t obey the third rule." % i)
        # print(dictionary)
    return dictionary


if __name__ == '__main__':
    """Call run_assembly_tester(). Do not modify this block."""
    run_assembly_tester()
