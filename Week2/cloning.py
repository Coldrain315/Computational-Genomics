
from compsci260lib import get_fasta_dict   # this will also import the sys and re modules


def solve_cloning():
    """Function demonstrating the procedure of extracting a gene and inserting
    into a plasmid using restriction enzymes."""

    aim2_fasta_path = 'aim2_plus_minus_1kb.fasta'
    pRS304_fasta_path = 'pRS304.fasta'



    # Read in the aim2 genomic sequence from the fasta file along with its
    # upstream and downstream regions.
    #
    # Your code goes here
    #
    aim2_genomic = ""

    # Store the beginning and the end of the Aim2 gene as Python indices.
    aim2_beg = None  # Your value goes here
    aim2_end = None  # Your value goes here

    # Report start, end, length in nucleotides, and length in amino acids
    #
    # Your code goes here
    #

    # Define regular expression terms for each restriction enzyme
    r_enzymes = get_restriction_enzymes_regex()

    # Store coordinates of restriction sites found upstream, downstream, and
    # within the aim2 gene
    r_enzyme_sites = find_aim2_restriction_enzymes(
        aim2_beg, aim2_end, aim2_genomic)

    # Report the found restriction enzyme sites
    #
    # Your code goes here
    #

    # Read in the pRS304 plasmid sequence and find restriciton sites
    # in the entire plasmid
    prs304_genomic = None # Your code here
    p_enzyme_sites = find_pRS304_restriction_sites(prs304_genomic)

    # Report relevant summary information
    mcs_start = None # Your value here
    mcs_end = None # Your value here
    mcs_enzyme_sites = report_pRS304_MCS_sites(p_enzyme_sites, mcs_start,
        mcs_end)

    # Extract aim2 gene and insert into the plasmid, report its length
    #
    # Your code goes here
    #


def get_restriction_enzymes_regex():
    """Returns a dictionary of restriction enzyme regular expressions for
    searching in genomic sequences.

    This function should be used for find_aim2_restriction_enzymes and
    find_pRS304_MCS_restriction_sites.
    """

    r_enzymes = {
        "BamHI": "<regular expression to identify BamHI site>",
        #
        # ...
        #
    }
    return r_enzymes


def find_aim2_restriction_enzymes(aim2_beg, aim2_end, aim2_genomic):
    """Finds the restriction enzyme sites in the aim2 genomic sequence. Stored
    as a dictionary of lists. Each restriction enzyme key corresponds to a list
    of dictionaries containing relevant information for a found site.

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched in the genome
        start (int): the start nucleotide position in the genome
        end (int): the ending position
        location (str): the position of the site relative the aim2 gene.
            (must be "upstream", "downstream", or "within")

    A valid returned dictionary may look like this:
    {
        'BamHI' : [
            {
                'start': 10,
                'end': 15,
                'sequence': 'GGATCC',
                'location': 'upstream'
            },
            {
                'start': 100,
                'end': 105,
                'sequence': 'GGATCC',
                'location': 'downstream'
            }
        ],
        'BstYI' : [
            {
                'start': 30,
                'end': 35,
                'sequence': 'AGATCC',
                'location': 'within'
            }
        ]
    }
    """

    # load restriction enzyme regular expressions
    r_enzymes = get_restriction_enzymes_regex()

    r_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}

    #
    # Your code goes here
    #

    return r_enzyme_sites



def report_pRS304_MCS_sites(p_enzyme_sites, mcs_start, mcs_end):
    """
    For each restriction enzyme, 
    
    Report how often that enzyme cuts the plasmid outside the MCS, how often
    it cuts the plasmid inside the MCS, and relevant details about any sites
    located inside the MCS.
    """
    
    #
    # Your code goes here
    #

    return # placeholder code, for a valid function


def find_pRS304_restriction_sites(prs304_genomic):
    """Finds the restriction sites in pRS304 genomic sequence. Stored as a
    dictionary of lists. Each restriction enzyme key corresponds to a list of
    dictionaries containing relevant information for a found site.

    This code will be similar (and simpler) code 
    found in `find_aim2_restriction_enzymes`. So it may be helpful to use code
    you wrote there.

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched
        start (int): the start nucleotide position in the plasmid
        end (int): the ending position in the plasmid

    A valid returned dictionary may look like this:
    {
        'BamHI' : [
            {
                'start': 10,
                'end': 15,
                'sequence': 'GGATCC'
            },
            {
                'start': 100,
                'end': 105,
                'sequence': 'GGATCC'
            }
        ],
        'BstYI' : [
            {
                'start': 30,
                'end': 35,
                'sequence': 'AGATCC'
            }
        ]
    }
    """

    # load restriction enzyme regular expressions
    r_enzymes = get_restriction_enzymes_regex()

    p_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}
    #
    # Your code goes here
    #

    return p_enzyme_sites


if __name__ == '__main__':
    """Run solve_cloning(). Do not modify this code"""
    solve_cloning()
