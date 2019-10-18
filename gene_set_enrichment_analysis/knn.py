import sys
import pandas as pd 
import math

def KNN(expression_df, samples_df, k, p):
	"""
		Runs the k-nearest neighbor algorithm on gene expression data to predict whether a patient is sick or healthy

		Inputs:
		   expression_df = pandas data frame containing gene expression for each sample
		   samples_df = pandas data frame containing samples and their classification (1 = patient, 0 = healthy)
		   k = number of neighbors to consider for k-nearest neighbor
		   p = the fraction of neighbors needed for a positive classification
		Returns:
		   combined_df = a pandas data frame that contains the sample name and the predicted classification (1 = patient, 0 = healthy) for the sample
		   sensitivity = true positives / (true positives + true negatives) when comparing knn predictions to ground truths
		   specificity = true negatives / (true negatives + false positives) when comparing knn predictions to ground truths
	"""
	classifications = pd.DataFrame(index = range(0,len(samples_df.index)), columns = range(0,len(samples_df.columns)))
	classifications=classifications.fillna(0)

	for col in range(1,len(expression_df.columns)):
		nearest_neigh = {}
		neighbors = []

		classifications.iloc[[col-1],[0]] = expression_df.columns[col]
		
		# Calculates euclidean distance between each subject
		for j in range(1,len(expression_df.columns)):
			distance = 0
			for i in range(0,len(expression_df.index)):
				distance += (expression_df.iat[i,col] - expression_df.iat[i,j])**2
			euc_dist = math.sqrt(distance)
			if euc_dist != 0:
				nearest_neigh[expression_df.columns[j]]=euc_dist

		# Finds the k nearest neighbors (closest euclidean distance)
		for x in range(k): 
			key = min(nearest_neigh, key=nearest_neigh.get)
			neighbors.append(key)
			nearest_neigh.pop(key, None)
		
		# Makes prediction of healthy or sick depending on if enough neighbors are sick patients (exceed threshold)
		# (sum outcomes / number of neighbors) > threshold
		sum_outcomes = 0
		for neighbor in neighbors:
			sum_outcomes += samples_df[samples_df[0] == neighbor].iat[0,1]

		if (sum_outcomes/k) > p:
			classifications.iloc[[col-1],[1]] = 1
		else:
			classifications.iloc[[col-1],[1]] = 0

	combined_df = samples_df.join(classifications.set_index(0), on=0, lsuffix='Sample', rsuffix='Prediction')
	sensitivity, specificity = calc_sens_spec(combined_df)

	classifications = classifications.sort_values(by=[0])

	return classifications, sensitivity, specificity

def calc_sens_spec(combined_df):
	"""
		Creates a confusion matrix by counting the number of tp, fp, tn, and fn. 
		Then, returns sensitivity and specificity

		Inputs:
		   combined_df = a pandas data frame that contains three columns: sample name, sample classification, KNN predictions
		Returns:
		   sensitivity = true positives / (true positives + true negatives)
		   specificity = true negatives / (true negatives + false positives)
	"""
	tp = 0	# true positive
	fp = 0	# false positive
	tn = 0	# true negative
	fn = 0	# false negative
	for row in range(len(combined_df.index)):
		if combined_df['1Sample'][row] == 1:
			if combined_df['1Sample'][row] == combined_df['1Prediction'][row]:
				tp += 1
			else:
				fn += 1
		else:
			if combined_df['1Sample'][row] == combined_df['1Prediction'][row]:
				tn += 1
			else:
				fp += 1
	
	# returns: sensitivity, specificity
	return tp/(tp+fn), tn/(tn+fp)

def main():

	# Check that the file is being properly used
	if (len(sys.argv) != 5):
		print("Please specify an expression file, samples file, k, and p as args.")
		return
		
	# Input files
	expression_file = sys.argv[1]
	samples_file = sys.argv[2]
	k = int(sys.argv[3])
	p = float(sys.argv[4])

	# Load files into pandas data frames
	expression_df = pd.read_csv(expression_file, sep='\t')
	samples_df = pd.read_csv(samples_file, sep='\t', header=None)

	# Runs KNN 
	classifications, sensitivity, specificity = KNN(expression_df, samples_df, k, p)

	# Outputs the results from KNN into a txt file
	with open('sample_assignments.txt', 'w') as output_f:
		for row in range(len(classifications.index)):
			output_f.write(classifications.iat[row,0] + '\t' + str(classifications.iat[row,1]) + '\n')

	# Prints the sensitivity and specificity
	print()
	print('Sensitivity: ', "%.2f" % round(sensitivity,2), '\n')
	print('Specificity: ', "%.2f" % round(specificity,2))


if __name__=="__main__":
	main()