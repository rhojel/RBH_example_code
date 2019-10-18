import sys
import numpy as np
from math import exp, ceil, log 
import csv

def hashFunc(a, b, n_buckets, x, p = 123457):
	y = x % p
	hash_val = (a*y+b) % p
	return int(hash_val % n_buckets)


def buildHashTable(delta, n_buckets, hashParam, word_stream):

	F = np.zeros((int(log(1/delta)),n_buckets))
	t = 0

	with open(word_stream, 'r') as in_f:
		reader = csv.reader(in_f, delimiter = '\t')
		for row in reader:
			t += 1
			for j in range(int(log(1/delta))):
				loc = hashFunc(hashParam[j,0],hashParam[j,1], n_buckets, int(row[0]))
				F[j,loc] += 1

	return F, t

def calcError(delta, hashParam, counts_file, approx_F, n_buckets):

	error = []
	actual = []
	counts = np.zeros(int(log(1/delta)))

	with open(counts_file, 'r') as in_f:
		reader = csv.reader(in_f, delimiter = '\t')
		for row in reader:
			for j in range(int(log(1/delta))):
				loc = hashFunc(hashParam[j,0],hashParam[j,1], n_buckets, int(row[0]))
				counts[j] = approx_F[j,loc]

			approx_count = np.amin(counts)

			error.append((approx_count - int(row[1]))/int(row[1]))
			actual.append(int(row[1]))

	return error, actual

def main():

	delta = exp(-5)
	epi = exp(1) * 10**-4

	n_buckets = ceil(exp(1)/epi)
	word_stream = sys.argv[1]
	counts_file = sys.argv[2]
	hash_params = sys.argv[3]

	hashParam = np.zeros((int(log(1/delta)),2))

	with open(hash_params, 'r') as in_f:
		reader = csv.reader(in_f, delimiter = '\t')
		i = 0
		for row in reader:
			hashParam[i,] = (int(row[0]),int(row[1]))
			i += 1

	approx_F, t = buildHashTable(delta, n_buckets, hashParam, word_stream)

	error, actual = calcError(delta, hashParam, counts_file, approx_F, n_buckets)

	actual = np.asarray(actual)/t

	with open("Ouput_Error_Full.txt", 'w') as o_f:
		for i in range(len(error)):
			o_f.write(str(error[i])+"\t"+ str(actual[i])+"\n")

	# with open("Ouput_Actual_2", 'w') as o_f:
	# 	for i in range(len(actual)):
	# 		o_f.write(str(actual[i])+",")


if __name__=="__main__":
	main()

## Graphing done in Jupyter Notebook
## Uncomment to Run
# import csv
# import numpy as np
# import matplotlib.pyplot as plt

# error = []
# actual = []

# with open("Ouput_Error_Full.txt", 'r') as in_f:
#     reader = csv.reader(in_f, delimiter = '\t')
#     for row in reader:
#         error.append(float(row[0]))
#         actual.append(float(row[1]))

# error = np.asarray(error)
# actual = np.asarray(actual)

# plt.loglog(actual,error,'.',markersize=2)
# plt.xlabel("Log Exact Frequencies")
# plt.ylabel("Log Relative Error")
# plt.title("Exact Word Frequency vs Relative Count Error \n \
#           for Data Stream Count Approximation Algorithm \n \
#           # of hash func = 5, n_buckets = 10^4")
