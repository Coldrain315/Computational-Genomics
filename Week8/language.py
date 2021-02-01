import textwrap, random
from compsci260lib import *


def run_language():
    """Call `solve_language` with the provided parameters."""

    file_name1 = "tidy.heart.of.darkness.txt"
    file_name2 = "tidy.paradise.lost.txt"
    num_outtext = 1
    for i in range(4):
        order = i+1
        solve_language(file_name1, order, num_outtext)
        solve_language(file_name2, order, num_outtext)


def solve_language(file_name, order, num_outtext):
    """Create a Markov model for a given text file and output artificially
    generated text from the model.

    Args:
        file_name (str): path of the text to process
        order (int): order of the Markov model
        num_outtext (int): number of output texts to generate
    """

    # Read the contents of the file
    f = open(file_name, 'r', encoding="utf-8")

    if f is None:
        print("Can't open " + file_name)
    else:
        contents = f.read()
        f.close()
        contents = contents.replace("\n", "")
        contents = contents.replace("\r", "")
    # print(contents)
    # This dictionary will store all the data needed to estimate the Markov model:
    # contents = "agggcagcgggcg"
    # Step 1: Count up occurrences of the various k-tuples in the training text
    print("Building dict of k-tuples and their counts...")
    # txt_dict = build_dict(contents, order)
    txt_dict = build_dict(contents, order)
    display_dict(txt_dict)
    # Step 2: Collect the counts necessary to estimate transition probabilities
    print("Collecting the counts to estimate transition probabilities...")
    txt_dict = collect_counts(contents, order, txt_dict)
    display_dict(txt_dict)

    # Step 3: Generate artificial text from the trained model
    for _ in range(num_outtext):
        seed = contents[0:order]
        M = 1000
        text = seed

        print("\nOne version of the story is:")

        # Display my story
        for _ in range(M):
            next_character = generate_next_character(seed, txt_dict)
            text += next_character
            seed = seed[1:] + next_character

        text_list = textwrap.wrap(text, 72)
        text = "\n".join(text_list)
        print(text)

    # return the last generated
    return text


def display_dict(txt_dict):

    print("key\tcount\tfollowers")
    for key in sorted(txt_dict):

        # k-tuple to int
        if type(txt_dict[key]) == int:
            print("%s\t%d "%(key, txt_dict[key]) , "\t", end=' ')
        # k-tuple to dict
        else:
            print("{} {:>6} ".format(key, txt_dict[key]["count"]), end='\t')
            for f, fcount in sorted(txt_dict[key]["followers"].items()):  # let the displayed dictionary be sorted
                print('%s:%d' % (f, fcount), end=' ')
        print()


def build_dict(contents, k):
    """Builds a dictionary of k-character (k-tuple) substring counts. Store the
    dictionary mapping from the k-tuple to an integer count.

    Args:
        contents (str): the string contents of to count
        k (int): number of characters in the substring

    Returns:
        a text dictionary mapping k-tuple to an integer. Example output with
        k=2:
        { 
            'ac': 1, 
            'cg': 2, 
            ... 
        }
    """

    #
    dictk = {}
    for i in range(len(contents) - 1):
        if i + k < len(contents):
            # the current substring; the key
            substr = contents[i:i + k]
        else:
            break
        if substr in dictk.keys():  # if the first time the key appears, initialize the key and its value
            dictk[substr] += 1
        else:
            dictk[substr] = 1  # if the key already exists, add 1 to its count
    return dictk


def collect_counts(contents, k, txt_dict):
    """Modify k-tuple dictionary to a mapping from k-tuple to a dictionary of
    of counts and dictionary of follower counts.
    
    Args:
        contents (str): the string contents of to count
        k (int): number of characters in the substring
        txt_dict (dict): mutable dictionary to store the k-character counts

    Returns;
        a dictionary mapping k-tuple to a dictionary of counts and dictionary
        of follower counts. Example:
        {
            'ac': {
                'count': 1,
                'followers': {'g': 1, 'c': 2}
            },
            ...
        }
        
    Note: This function is similar and an extension of `build_dict`. We 
    separated the k-character and follower counting to explain each as 
    distinct concepts. It is possible to implement `collect_counts` without
    using the txt_dict returned from `build_dict`.

    Our auto-grader will test the behavior of `build_dict` so that function 
    does need to work properly, but when you actually go to train your model 
    on the provided texts, you can simply call collect_counts() to save the
    time. 

    """

    #
    dictk = {}
    for i in range(len(contents) - 1):
        if i + k < len(contents):
            substr = contents[i:i + k]  # the current substring; the key at this time
            nextstr = contents[i + k]  # the subsequent char of current string
        else:
            break
        if substr in dictk.keys():   # if the current substring already exist as key, add 1 to its count
            dictk[substr]['count'] += 1
            if nextstr in dictk[substr]['followers'].keys():  # if the next string already exists as key, add 1 to its count
                dictk[substr]['followers'][nextstr] += 1
            elif nextstr not in dictk[substr]['followers'].keys():  # if the first time the next char of current string appears, initialize the key and its value
                dictk[substr]['followers'][nextstr] = 1
        else:
            dictk[substr] = {}  # if the first time the key appears, initialize the key and its values, including the count and its follower
            dictk[substr]['count'] = 1
            dictk[substr]['followers'] = {}
            dictk[substr]['followers'][nextstr] = 1
    return dictk


def generate_next_character(seed, txt_dict):
    """Randomly select the next character of a k-tuple using the follower
    counts to determine the probability.

    Args:
        seed (str): k-tuple to follow from
        txt_dict (dict): k-tuple count follower dictionary

    Returns:
        (str) of the next character
    """

    #
    nextstr = txt_dict[seed]
    rate_list = list(nextstr['followers'].values())
    returnstr = list(nextstr['followers'].keys())[random_index(rate_list)]
    return returnstr


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


if __name__ == '__main__':
    """Main method call, do not modify"""
    run_language()

