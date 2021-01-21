from compsci260lib import *  # this will also import the sys and re modules
import re


def solve_plasmid():
    # Load the plasmid reads from fasta ensuring they are
    # in the expected order to be assembled
    #
    # Your code here
    #
    plasmid_path = 'plasmid.fasta'
    plasmid_genome = get_fasta_dict(plasmid_path)
    # print(plasmid_genome)
    reads = [plasmid_genome['start']]
    for i in plasmid_genome:
        if i == 'start':
            continue
        else:
            reads.append(plasmid_genome[i])
    # Assemble the reads into a single-stranded linearized version of the
    # plasmid:
    plasmid_dna = simple_assembler(reads)
    print("The assembled plasmid sequence is %s" % plasmid_dna)
    # Report the length of the assembled plasmid
    plasmid_length = len(plasmid_dna)
    print("There are %d nucleotides in the assembled sequence." % plasmid_length)

    # Search for ORFs in the reconstructed plasmid DNA and report your results
    negative_strand = nucleotide_match(plasmid_dna)
    orfs = find_orfs_circular_double_stranded(plasmid_dna, min_length_aa=50)
    # print(orfs)
    orf_length = []
    for i in orfs:
        orf_length.append(i['nlength'])
        print("The orf on %s strand with frame %d starts at %d and ends at %d, with stop codon %s. The length of nucleotide is %d and length of amino acid is %d. " % (i['strand'], i['frame'], i['start'], i['stop'], i['stopcodon'], i['nlength'], i['aalength']))
    print("There are %d orfs found." % len(orfs))
    max_index = orf_length.index(max(orf_length))
    start = orfs[max_index]['start']
    end = orfs[max_index]['stop']
    if orfs[max_index]['strand'] == '+':
        if start < end:
            largest_orf = plasmid_dna[start-1:end]
        else:
            largest_orf = plasmid_dna[start-1:]+plasmid_dna[0:end]
    elif orfs[max_index]['strand'] == '-':
        if start < end:
            largest_orf = negative_strand[start-1:end]
        else:
            largest_orf = negative_strand[start-1:]+plasmid_dna[0:end]
    amino_acid = translate(largest_orf)
    print("The longest orf found has %d nucleotides" % max(orf_length))
    print("The translated amino acid is %s" % amino_acid)
    fraction = frac_orf(plasmid_dna, orfs)  # + frac_orf(negative_strand, orfs)
    print("The fraction of the genome that is coding is approximately %.2f" % fraction)
    # print_abbrv(plasmid_dna)


def simple_assembler(reads):
    """Given a list of reads, as described in the problem, return an assembled
    DNA sequence, as a string. For consistency, use the first entry in the
    fragment reads list as the starting position of the returned sequence.

    For example, if we were to take in a list of three reads, 31 nucleotides
    long each. The last 15 nucleotides of each read would overlap with one
    other read, and the assembled sequence would be 48 nucleotides long with
    the sequence starting with the beginning of the first read.

    Args:
        reads (list): list of sequence strings as reads

    Returns:
         str: an assembled genomic sequence as a string starting with the first
              read in `reads`
    """
    #

    i = reads[0]  # select the first entry as the beginning target read
    new_sequence = i
    temp_end = i[-15:]  # slice out the last 15 nucleotides for comparison
    new_reads = reads  # a new list that stores the reads, and when one read matches its start and end with other two reads, move the read out of the list
    while True:
        if len(
                new_reads) == 0:  # when there are no reads in the new list, it means all the reads have found their previous one and next one
            break
        else:
            for j in new_reads:
                temp_start = j[0:15]
                if temp_end == temp_start and len(
                        new_reads) > 1:  # search in the new read list, find the read with the first 15 nucleotides of a read matches the end of our target read
                    new_part = j[
                               15:]  # the nucleotide sequence of the found read, except for the beginning 15 nucleotides, because they already exist in the sequence as the end of previous read.
                    new_sequence += new_part  # add the found read to our sequence
                    new_reads.remove(j)  # remove the target read from the list since it already finds its match
                    i = j  # define the found one that matches the end of our target read as new target read
                    temp_end = i[-15:]
                    break
                elif temp_end == temp_start and len(new_reads) <= 1:
                    new_reads.remove(j)
                    new_sequence = new_sequence[0:-15]
                    break
    return new_sequence


def find_orfs_circular_double_stranded(seq, min_length_aa=0):
    """This is a function for finding sufficiently long ORFs in all reading
    frames in a sequence of double-stranded circular DNA.

    The function takes as input parameters: a string representing a genomic
    sequence, the minimum length (in amino acids) for an ORF before it will be
    returned (which defaults to 0).

    Args:
        seq (str): a genomic sequence
        min_length_aa (int): minimum length of found ORFs in amino acids

    Returns:
        list: of dictionaries with information on each ORF found.

    Where each ORF found is represented by a dictionary with
    the following keys:
        frame (int): The nucleotide offset in which the ORF was found. (Must be
        0, 1, or 2)
        start (int): the nucleotide position of the start of the ORF 
            (relative to 5' end of the found ORF's strand)
        stop (int): the nucleotide position of the end of the ORF
            (relative to 5' end of the found ORF's strand)
        stopcodon (str): the nucleotide triplet of the stop codon
        nlength (int): the length (in nucleotides) of the ORF
        strand (str): the strand of the found ORF (Must be '+' or '-')

    A valid return list may look something like this:
    [
        {
            'frame': 0,
            'stop': 13413,
            'aalength': 4382,
            'start': 265,
            'stopcodon': 'UAA',
            'nlength': 13149,
            'strand': '+'
        },
        {
            'frame': 0,
            'stop': 27063,
            'aalength': 221,
            'start': 26398,
            'stopcodon': 'UAA',
            'nlength': 666,
            'strand': '-'
        }
    ]
    """

    #
    positive = seq
    reverse_negative = nucleotide_match(positive)  # find the reverse-negative strand
    min_length_nt = min_length_aa * 3  # we are given the min length of amino acid, and we can get the min length of nucleotide which is three times of amino acid
    orfs = []  # the list that will store all orfs

    for frame in range(3):  # loop through three potential
        orf_list = []
        start_regex = re.compile(r"ATG", re.IGNORECASE)
        stop_regex = re.compile(r"T(AA|AG|GA)", re.IGNORECASE)  # using regular expression to find start and stop codon

        # we first go through the positive strand
        start_index = frame
        end_index = start_index
        in_orf = False  # is false when we haven't started or found a start codon
        orf_dict = {'frame': frame, 'strand': '+'}
        while start_index <= len(positive) and end_index <= 2 * len(positive):
            if not in_orf:
                if start_regex.match(positive, start_index):  # if the nucleotide at current position matches the start codon
                    in_orf = True  # turn it to true when we started the orf
                    orf_dict['start'] = start_index + 1  # 0-index to 1-index
                    end_index = start_index + 3
                start_index += 3
            else:  # this must be in an ORF
                if end_index <= len(positive):  #
                    if stop_regex.match(positive, end_index):
                        in_orf = False
                        nuc_length = end_index + 3 + 1 - orf_dict['start']
                        if nuc_length - 3 >= min_length_nt:
                            orf_dict['stop'] = end_index + 3
                            orf_dict['nlength'] = nuc_length
                            orf_dict['aalength'] = (nuc_length - 3) // 3
                            orf_dict['stopcodon'] = positive[end_index:end_index + 3]
                            orf_list.append(dict(orf_dict))
                    if end_index + 3 <= len(positive):
                        start_index = end_index + 3
                    else:
                        start_index = end_index
                    # print(start_index)
                elif end_index > len(positive):  # if end_index is larger than length of sequence, it means we are
                    # print(end_index)
                    used_end_index = end_index - len(positive)
                    if stop_regex.match(positive, used_end_index):
                        start_index = end_index
                        in_orf = False
                        nuc_length = end_index + 3 + 1 - orf_dict['start']
                        if nuc_length - 3 >= min_length_nt:
                            orf_dict['stop'] = used_end_index + 3
                            orf_dict['nlength'] = nuc_length
                            orf_dict['aalength'] = (nuc_length - 3) // 3
                            orf_dict['stopcodon'] = positive[used_end_index:used_end_index + 3]
                            orf_list.append(dict(orf_dict))
                end_index += 3
            # print(orf_list)

        # we then go through the negative strand; the code should have exactly same logic with the positive one
        start_index = frame
        end_index = start_index
        in_orf = False
        orf_dict = {'frame': frame, 'strand': '-'}
        # negative strand
        while start_index <= len(reverse_negative) and end_index <= 2 * len(reverse_negative):
            if not in_orf:
                if start_regex.match(reverse_negative, start_index):
                    in_orf = True
                    orf_dict['start'] = start_index + 1
                    end_index = start_index + 3
                start_index += 3
            else:  # we must be in an ORF
                if end_index <= len(reverse_negative):
                    if stop_regex.match(reverse_negative, end_index):
                        in_orf = False
                        nuc_length = end_index + 3 + 1 - orf_dict['start']
                        if nuc_length - 3 >= min_length_nt:
                            orf_dict['stop'] = end_index + 3
                            orf_dict['nlength'] = nuc_length
                            orf_dict['aalength'] = (nuc_length - 3) // 3
                            orf_dict['stopcodon'] = reverse_negative[end_index:end_index + 3]
                            orf_list.append(dict(orf_dict))
                    if end_index + 3 <= len(reverse_negative):
                        start_index = end_index + 3
                    else:
                        start_index = end_index
                    # print(start_index)
                elif end_index > len(reverse_negative):
                    # print(end_index)
                    used_end_index = end_index - len(reverse_negative)
                    if stop_regex.match(reverse_negative, used_end_index):
                        start_index = end_index
                        in_orf = False
                        nuc_length = end_index + 3 + 1 - orf_dict['start']
                        if nuc_length - 3 >= min_length_nt:
                            orf_dict['stop'] = used_end_index + 3
                            orf_dict['nlength'] = nuc_length
                            orf_dict['aalength'] = (nuc_length - 3) // 3
                            orf_dict['stopcodon'] = reverse_negative[used_end_index:used_end_index + 3]
                            orf_list.append(dict(orf_dict))
                end_index += 3
        orfs.extend(orf_list)
    return orfs


def nucleotide_match(seq):  # conduct wrack crick process
    positive = seq
    negative = ''
    for i in positive:
        if i == 'a':
            negative += 't'
        elif i == 't':
            negative += 'a'
        elif i == 'g':
            negative += 'c'
        elif i == 'c':
            negative += 'g'
    reverse_negative = negative[::-1]
    return reverse_negative


def frac_orf(seq, orfs_dict):
    len_seq = len(seq)
    lst_seq = [0 for i in range(len_seq)]
    lst_seq_wc = [0 for i in range(len_seq)]
    for orf in orfs_dict:
        if orf["stop"] < orf["start"]:
            if orf["strand"] == "+":
                for i in range(orf["stop"]-3):
                    lst_seq[i] = 1
                for i in range(orf['start']-1, len_seq):
                    lst_seq[i] = 1
            else:
                for i in range(orf["stop"]-3):
                    lst_seq_wc[i] = 1
                for i in range(orf['start'] - 1, len_seq):
                    lst_seq_wc[i] = 1
        else:
            if orf["strand"] == "+":
                for i in range(orf["start"] - 1, orf["stop"]-3):
                    lst_seq[i] = 1
            else:
                for i in range(orf["start"] - 1,orf["stop"]-3):
                    lst_seq_wc[i] = 1
                lst_seq_wc[orf["start"]-1:orf["stop"]] = [1 for i in range(orf["stop"]-orf["start"]+1)]
    return (sum(lst_seq)+sum(lst_seq_wc))/(len_seq*2)


if __name__ == '__main__':
    solve_plasmid()
