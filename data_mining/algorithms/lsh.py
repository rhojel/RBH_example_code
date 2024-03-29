# Authors: Jessica Su, Wanzi Zhou, Pratyaksh Sharma, Dylan Liu, Ansh Shukla

# Edited/Completed by Richard Hojel

import numpy as np
import random
import time
import pdb
import unittest
from PIL import Image
import sys
import matplotlib.pyplot as plt

# Finds the L1 distance between two vectors
# u and v are 1-dimensional np.array objects
def l1(u, v):
    return np.sum(np.abs(u-v))

# Loads the data into a np array, where each row corresponds to
# an image patch -- this step is sort of slow.
# Each row in the data is an image, and there are 400 columns.
def load_data(filename):
    return np.genfromtxt(filename, delimiter=',')

# Creates a hash function from a list of dimensions and thresholds.
def create_function(dimensions, thresholds):
    def f(v):
        boolarray = [v[dimensions[i]] >= thresholds[i] for i in range(len(dimensions))]
        return "".join(map(str, map(int, boolarray)))
    return f

# Creates the LSH functions (functions that compute L K-bit hash keys).
# Each function selects k dimensions (i.e. column indices of the image matrix)
# at random, and then chooses a random threshold for each dimension, between 0 and
# 255.  For any image, if its value on a given dimension is greater than or equal to
# the randomly chosen threshold, we set that bit to 1.  Each hash function returns
# a length-k bit string of the form "0101010001101001...", and the L hash functions 
# will produce L such bit strings for each image.
def create_functions(k, L, num_dimensions=400, min_threshold=0, max_threshold=255):
    functions = []
    for i in range(L):
        dimensions = np.random.randint(low = 0, 
                                   high = num_dimensions,
                                   size = k)
        thresholds = np.random.randint(low = min_threshold, 
                                   high = max_threshold + 1, 
                                   size = k)

        functions.append(create_function(dimensions, thresholds))
    return functions

# Hashes an individual vector (i.e. image).  This produces an array with L
# entries, where each entry is a string of k bits.
def hash_vector(functions, v):
    return np.array([f(v) for f in functions])

# Hashes the data in A, where each row is a datapoint, using the L
# functions in "functions."
def hash_data(functions, A):
    return np.array(list(map(lambda v: hash_vector(functions, v), A)))

# Retrieve all of the points that hash to one of the same buckets 
# as the query point.  Do not do any random sampling (unlike what the first
# part of this problem prescribes).
# Don't retrieve a point if it is the same point as the query point.
def get_candidates(hashed_A, hashed_point, query_index):
    return filter(lambda i: i != query_index and \
        any(hashed_point == hashed_A[i]), range(len(hashed_A)))

# Sets up the LSH.  You should try to call this function as few times as 
# possible, since it is expensive.
# A: The dataset.
# Return the LSH functions and hashed data structure.
def lsh_setup(A, k = 24, L = 10):
    functions = create_functions(k = k, L = L)
    hashed_A = hash_data(functions, A)
    return (functions, hashed_A)

# Run the entire LSH algorithm
def lsh_search(A, hashed_A, functions, query_index, num_neighbors = 10, return_dist = False):
    hashed_point = hash_vector(functions, A[query_index, :])
    candidate_row_nums = get_candidates(hashed_A, hashed_point, query_index)
    
    distances = map(lambda r: (r, l1(A[r], A[query_index])), candidate_row_nums)
    best_neighbors = sorted(distances, key=lambda t: t[1])[:num_neighbors]

    if return_dist:
        return [t[1] for t in best_neighbors]

    return [t[0] for t in best_neighbors]

# Plots images at the specified rows and saves them each to files.
def plot(A, row_nums, base_filename):
    for row_num in row_nums:
        patch = np.reshape(A[row_num, :], [20, 20])
        im = Image.fromarray(patch)
        if im.mode != 'RGB':
            im = im.convert('RGB')
        im.save(base_filename + "-" + str(row_num) + ".png")

# Finds the nearest neighbors to a given vector, using linear search.
def linear_search(A, query_index, num_neighbors, return_dist = False):
    distances = map(lambda r: (r, l1(A[r], A[query_index])), [i for i in range(np.size(A,0)) if i != query_index])
    best_neighbors = sorted(distances, key=lambda t: t[1])[:num_neighbors]

    if return_dist:
        return [t[1] for t in best_neighbors]

    return [t[0] for t in best_neighbors]

# A function that computes the error measure
def error_measure(A, k, L):
    functions, hash_A = lsh_setup(A, k, L)
    error = 0

    for index in range(99,1000,100):
        dist_lsh = []
        while len(dist_lsh) != 3:
            dist_lsh = lsh_search(A, hash_A, functions, index, 3, True)
        
        dist_linear = linear_search(A, index, 3, True)
        error += (np.sum(dist_lsh)/np.sum(dist_linear))

    return error/10


# Solve Problem 4
def problem4():
    # Load data
    A = load_data(sys.argv[1])
    # Create hash and hash functions
    functions, hash_A = lsh_setup(A)

    # Calculates the avg runtime of finding the 3 nearest neighbors for both algorithms
    # Part i (Comment out if not running this part)
    lsh_total = 0
    linear_total = 0
    for index in range(99,1000,100):
        three_neighbors_lsh = []
        while len(three_neighbors_lsh) != 3:
            lsh_start = time.time()
            three_neighbors_lsh = lsh_search(A, hash_A, functions, index, 3)
            lsh_end = time.time()

        linear_start = time.time()
        three_neighbors = linear_search(A, index, 3)
        linear_end = time.time()

        lsh_total += (lsh_end - lsh_start)
        linear_total += (linear_end - linear_start)

    lsh_avg = lsh_total / 10
    linear_avg = linear_total / 10
    print('Average Search Time of LSH (3 neighbors): ' + str(lsh_avg)+'\n')
    print('Average Search Time of Linear Search (3 neighbors): ' + str(linear_avg)+'\n')

    # Calculating the error measure
    # Part ii (Comment out if not running this part)
    L_error = []
    k_error = []

    for L in range(10,22,2):
        L_error.append(error_measure(A,24,L))

    for k in range(16,26,2):
        k_error.append(error_measure(A,k,10))

    # Plot of L error
    plt.plot(L_error)
    plt.title('Error Value as a Function of L')
    plt.ylabel('Error Value')
    plt.xlabel('Values of L')
    locs, labels = plt.xticks()
    plt.xticks(locs,range(10,22,2))
    plt.savefig('L_error.png')
    plt.clf()

    # Plot of k error
    plt.plot(k_error)
    plt.title('Error Value as a Function of k')
    plt.ylabel('Error Value')
    plt.xlabel('Values of k')
    locs, labels = plt.xticks()
    plt.xticks(locs,range(16,26,2))
    plt.savefig('k_error.png')

    # Plotting the top 10 neighbors for row 100
    # Part iii (Comment out if not running this part)
    neighbors_linear = linear_search(A, 99, 10)
    
    neighbors_lsh = []
    while len(neighbors_lsh) != 10:
        neighbors_lsh = lsh_search(A, hash_A, functions, 99, 10)

    plot(A, [99], 'Original_Image')
    plot(A, neighbors_linear, 'Linear_Search_Neighbor')
    plot(A, neighbors_lsh, 'LSH_Neighbor')

#### TESTS #####

class TestLSH(unittest.TestCase):
    def test_l1(self):
        u = np.array([1, 2, 3, 4])
        v = np.array([2, 3, 2, 3])
        self.assertEqual(l1(u, v), 4)

    def test_hash_data(self):
        f1 = lambda v: sum(v)
        f2 = lambda v: sum([x * x for x in v])
        A = np.array([[1, 2, 3], [4, 5, 6]])
        self.assertEqual(f1(A[0,:]), 6)
        self.assertEqual(f2(A[0,:]), 14)

        functions = [f1, f2]
        self.assertTrue(np.array_equal(hash_vector(functions, A[0, :]), np.array([6, 14])))
        self.assertTrue(np.array_equal(hash_data(functions, A), np.array([[6, 14], [15, 77]])))

    ### TODO: Write your tests here (they won't be graded, 
    ### but you may find them helpful)


if __name__ == '__main__':
#    unittest.main() ### TODO: Uncomment this to run tests
    problem4()
