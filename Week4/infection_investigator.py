from bwt import *
from read_aligner import *
from compsci260lib import *


def align_patient_reads():
    """Align patient reads to each bacterial genome"""

    # Patients
    patients = ['patient1',
                'patient2',
                'patient3']

    # Bacterial species
    panel = ['Bacteroides_ovatus',
             'Bacteroides_thetaiotaomicron',
             'Bifidobacterium_longum',
             'Eubacterium_rectale',
             'Lactobacillus_acidophilus',
             'Peptoniphilus_timonensis',
             'Prevotella_copri',
             'Roseburia_intestinalis',
             'Ruminococcus_bromii',
             'Vibrio_cholerae'
             ]

    # Store the genome sequence and bwt structures for each bacterial species
    bact_sequences = {}
    bact_bwt_structures = {}
    for bacteria in panel:
        # For each of the bacteria in our panel, create a dictionary entry
        # where keys are bacterial names and values are their genome sequences.
        bact_sequences[bacteria] = list(get_fasta_dict(
            'reference_genomes/%s.fasta' % bacteria).values())[0]

        # For each of the bacterial genome sequences, create a dictionary entry
        # where keys are bacterial names and values are their BWT data
        # structures.
        bact_bwt_structures[bacteria] = make_all(bact_sequences[bacteria])

    # Store a mapping of patient names to reads
    patient_reads = {}
    for patient in patients:
        reads = list(get_fasta_dict('patients/%s.fasta' % patient).values())
        patient_reads[patient] = reads

    # Store the reads mapped per bacterial species per patient
    mapped_patient_reads = {}

    # Consider all patients
    for patient in patients:
        # Find uniquely mapped reads for each bacteria for the patient
        # and store the read start positions
        mapped_reads = find_aligned_reads_for_patient(
            patient_reads[patient], bact_bwt_structures)
        # print(mapped_reads)
        mapped_patient_reads[patient] = mapped_reads

    # Report the microbe prevalences for each of your patients
    for patient in patients:
        patient_read = mapped_patient_reads[patient]
        panel_dict = {}
        count = 0
        for bacteria in panel:
            panel_dict[bacteria] = len(patient_read[bacteria])
            count += len(patient_read[bacteria])
        for new_bacteria in panel:
            panel_dict[new_bacteria] = panel_dict[new_bacteria] / count
            print("For the %s, the prevalence of bacteria %s is %d" % (patient, new_bacteria, panel_dict[new_bacteria]))

    # Use `read_mapper` and `longest_zeroes` to identify unmapped regions
    # for the relevant patients and species (questions 2f-h)
    #
    for patient in patients:
        patient_read = mapped_patient_reads[patient]
        for bacteria in panel:
            mappings = read_mapper(patient_read[bacteria], bact_sequences[bacteria])
            if longest_zeros(mappings) is not None:
                start_pos, end_pos = longest_zeros(mappings)
                print("For the %s, the longest_zero on bacteria %s is from %d to %d" % (
                    patient, bacteria, start_pos, end_pos))
            else:
                print("There was no run of zeroes found.")

    # The part here below reports specifically patient1 and patient2 on Vibrio_cholerae
    # patient_read1 = mapped_patient_reads["patient1"]
    # mappings1 = read_mapper(patient_read1["Vibrio_cholerae"], bact_sequences["Vibrio_cholerae"])
    # start_pos1, end_pos1 = longest_zeros(mappings1)
    # patient_read2 = mapped_patient_reads["patient2"]
    # mappings2 = read_mapper(patient_read2["Vibrio_cholerae"], bact_sequences["Vibrio_cholerae"])
    # start_pos2, end_pos2 = longest_zeros(mappings2)
    #
    # print("For the patient1, the longest_zero on Vibrio_cholerae is from %d to %d" % (start_pos1, end_pos1))
    # print("For the patient2, the longest_zero on Vibrio_cholerae is from %d to %d" % (start_pos2, end_pos2))


def find_aligned_reads_for_patient(reads, bact_bwt_structures):
    """
    Given a list of reads for a patient, identify the reads that uniquely map
    to each bacteria's genome using its relevant BWT data structure. Reads that
    are mapped to a bacteria's genome should be stored as starting positions in
    a list.

    Args:
        reads (list of str):
            mapping a patient's read names (str) to read sequences (str).

        bact_bwt_structures: dictionary mapping bacterial names (str) to 
            structures required for efficient exact string matching). Refer
            to the return tuple from `make_all` in read_aligner.py.

    Returns:
        a dictionary mapping bacterial names (str) to the start positions of
        reads mapped to that bacteria's genome. Note that the end positions are
        not needed as the reads are all 50 bases long. Start positions should
        be stored using 0-indexing.

        Example return structure:
        {
            'Bacteroides_ovatus': [0, 100, 200, ...],
            ...
        }
    """
    reverse_reads = [reverse_complement(i) for i in reads]
    counts = {}
    temp_count = {}
    panel = list(bact_bwt_structures.keys())
    list_read_index = []
    indices_list = []
    # go through all reads of a patient
    for read_index, read in enumerate(reverse_reads):
        # go through all bacteria to find match
        for bact_index, bacterium in enumerate(panel):
            bwt_structures = bact_bwt_structures[bacterium]
            # apply find() function to find matches between reads and genome
            all_match_index = find(read, bwt_structures)
            # if matches found, store the position of first match
            if len(all_match_index) >= 1:
                match_index = all_match_index[0]
                # list_read_index stores the reads that are already matched with some other bacteria
                # check if current read in the list
                if read_index not in list_read_index:  # if not, it is valid, append it to list_read_index
                    list_read_index.append(read_index)
                    # add the read to dictionary of current bacterium
                    if bacterium in temp_count.keys():
                        temp_count[bacterium][read_index] = match_index
                    else:
                        temp_count[bacterium] = {read_index: match_index}
                else:  # if so, the read is not valid, find the previous bacterium that has match with this read
                    # and remove it from the dictionary of that bacterium
                    for prev_bact in range(bact_index):
                        if panel[prev_bact] in temp_count.keys():
                            judge = read_index in temp_count[panel[prev_bact]].keys()
                            if judge:
                                temp_count[panel[prev_bact]].pop(read_index)
                    # we can directly break the loop of bacterium for current read once the read is found matched with two bacterium
                    break
    for bact in panel:
        if bact in temp_count.keys():
            indices_list.append(list(temp_count[bact].values()))
        else:
            indices_list.append([])
    counts = dict(zip(panel, indices_list))
    return counts


def read_mapper(starting_positions, genome):
    """
    Using the starting positions of reads that were aligned to a bacterial
    genome, construct a count vector (as a list) that counts how many reads
    were aligned to a position in the genome. The vector's size will be the 
    length of the genome. You may assume that each read is 50 base pairs long.

    Args:
        starting_positions (list of ints): 
            starting positions of reads that were aligned to a genome.

        genome (str): the genomic nucleotide sequence to which the reads 
                      were aligned.

    Returns:
        (list of ints): vector of aligned read counts to the genome.
        i.e. [c_1, c_2, ..., c_i, ..., c_n], where n=length of the genome
        and c_i = the count of aligned reads for the patient at genome
        position i.
    """
    # initialize the count_vector as all 0
    count_vector = [0 for i in range(len(genome))]
    # find the indices that have starts
    for i in starting_positions:
        # add 1 to the number of matches at certain location i since it's in start location array
        count_vector[i] = count_vector[i] + 1
    return count_vector


# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.


def reverse_complement(seq):
    """Returns the reverse complement of the input string."""
    comp_bases = {'A': 'T',
                  'C': 'G',
                  'G': 'C',
                  'T': 'A'}
    rev_seq = list(seq)
    rev_seq = rev_seq[::-1]
    rev_seq = [comp_bases[base] for base in rev_seq]
    return ''.join(rev_seq)


def longest_zeros(count_vector):
    """Given a count vector, return the start and stop position of the longest
    string of zeros in the vector.

    Args:
        count_vector (list of ints): vector of aligned read counts to the
        genome.  see: the return value for `read_mapper()`

    Returns:
        (tuple of (int, int)): the start and stop position of the longest
        string of zeros in the count_vector.
    """

    zero_nums = []
    genome_len = len(count_vector)
    for i in range(0, genome_len):
        if count_vector[i] == 0:
            zero_nums.append(i)

    counter = 1
    longest_run = {}
    longest_run[1] = []
    longest_run[1].append(zero_nums[0])
    for z in zero_nums[1:]:
        if (z - longest_run[counter][-1]) == 1:
            longest_run[counter].append(z)
        else:
            counter += 1
            longest_run[counter] = []
            longest_run[counter].append(z)

    for run in list(longest_run.keys()):
        if 0 in longest_run[run] or genome_len - 1 in longest_run[run]:
            del longest_run[run]

    longest = []
    for run in list(longest_run.values()):
        if len(run) > len(longest):
            longest = run

    # There was no run of zeroes found, return None
    if len(longest) == 0:
        return None

    start = longest[0]
    stop = longest[-1]

    # Return the start and end positions of the longest zeroes (0-indexed)
    return start, stop


if __name__ == '__main__':
    align_patient_reads()
