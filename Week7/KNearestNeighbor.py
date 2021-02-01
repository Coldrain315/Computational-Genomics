import sys
from math import sqrt


def run_k_nearest_neighbor():
    """Run K Nearest neighbors against the training and test sets"""
    input_patients = read_training_data("gene_expression_training_set.txt")
    new_patients = read_test_data("gene_expression_test_set.txt")
    solve_k_nearest(input_patients, new_patients, k = 5)


def solve_k_nearest(input_patients, new_patients, k):
    """
    Read in the input patients as training data to use to determine and report
    the prognosis of the 10 new patients.

    Args:
        input_patients (list of dicts): dictionaries of training patient
            data. see: `read_training_data`

        new_patients (list of dicts): dictionaries of test patient data.
            see: `read_test_data`
    """
    
    #
    # YOUR CODE HERE
    #


def read_training_data(file_name):
    """Read the training gene expression data from a text file. Note: the
    patients in the training data are classified as "R" (responsive to
    treatment) or "N" (non-responsive to treatment).  For example,
    input_patients[0]["class"] = the class of the first patient (R or N)
    input_patients[0]["expression"][0] = the expression of the first
    gene for the first patient.

    Returns:
        (list of dicts): list of patients as a class and expression data. The
        dictionary of each patient will be in the form of:
            'class' -> string with values strictly 'N' or 'R' for
            non-responsive or responsive to the treatment
            'expression' -> list of floats of gene expression values

        and look something like:
            {'class': 'N', 'expression': [9.049, 8.313, ..., 6.428700888]}
    """
        
    f = open(file_name, "r")
  
    if f is None:
        sys.stderr("File " + file_name + " does not exist") 
    else:
        input_patients = [] 
        lines = f.readlines()
       
        for line in lines:
            line = line.strip()
            data = line.split()      # check if you are splitting on "\t" which is what you want to split on
            class_name = data.pop(0)
            float_data = [float(datum) for datum in data]   #convert to float data
            patient = {"class": class_name, "expression": float_data}
            input_patients.append(patient)
   
        f.close()
        return input_patients


def read_test_data(file_name):
    """Read the test gene expression data from a text file. Note: the
    patients in the test data are not classified.

   Returns:
    (list of dicts): list of patients as a class and expression data. The
    dictionary of each patient will be in the form of:
        'class' -> string with only 'unknown' as its value
        'expression' -> list of floats of gene expression values

    and look something like:
        {'class': 'unknown', 'expression': [9.049, 8.313, ..., 6.428700888]}
    """
    
    f = open(file_name, "r")
  
    if f is None:
        sys.stderr("File " + file_name + " does not exist") 
    else:
        test_patients = [] 
        lines = f.readlines()
       
        for line in lines:
            line = line.strip()
            data = line.split()      # check if you are splitting on "\t" which is what you want to split on
            float_data = [float(datum) for datum in data]
            patient = {"class": "unknown", "expression": float_data}
            test_patients.append(patient)
   
        f.close()
        return test_patients
     

def compute_dist(tuple_1, tuple_2):
    """Return the Euclidean distance between two points in any number of
    dimensions."""
    
    if len(tuple_1) != len(tuple_2):
        sys.stderr("Cannot compute Euclidean distance between tuples of different sizes!")
    
    dist = 0
    for i in range(len(tuple_1)):
        dist += (tuple_1[i] - tuple_2[i]) * (tuple_1[i] - tuple_2[i])
      
    return sqrt(dist)


if __name__ == '__main__':
    run_k_nearest_neighbor()
