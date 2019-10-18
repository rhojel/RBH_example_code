import sys
import pandas as pd 
import numpy as np
import math

def GSEA(expression_df, samples_df, gene_set_dict):
	"""
		Runs a gene set enrichment analysis and returns the enrichment score for each gene set

		Inputs:
		   expression_df = pandas data frame containing gene expression for each sample
		   samples_df = pandas data frame containing samples and their classification (1 = patient, 0 = healthy)
		   gene_set_dict = a dictionary mapping gene_set names to a set of genes in that gene_set
		Returns:
		   ES = a list of tuples tha contain the gene set name and the enrichment score of that gene set
	"""
	expression_df = expression_df.set_index('SYMBOL')
	healthy = samples_df[samples_df[1] == 0]
	sick = samples_df[samples_df[1] == 1]

	# Calculates the mean gene expression for sick patients and for healthy patients 
	mean_healthy = np.mean(expression_df[list(healthy[0])].values, axis =1)
	mean_sick = np.mean(expression_df[list(sick[0])].values, axis =1)
	gene_names = list(expression_df.index)

	# Gets the ratio of the sick patient expression to healthy patient expression and ranks the list from greatest ratio 
	# to smallest ratio
	ranked_genes = []
	for gene in range(len(gene_names)):
	    ranked_genes.append((gene_names[gene],mean_sick[gene]-mean_healthy[gene]))

	ranked_genes.sort(key=lambda tup: tup[1], reverse=True)

	# Conducts a random walk (specifically a brownian bridge) to find the enrichment score for each gene set.
	# Stores the gene set names and enrichment scores in a list
	ES = []
	for gene_set in gene_set_dict:
		up_score = math.sqrt((len(ranked_genes) - len(gene_set_dict[gene_set]))/len(gene_set_dict[gene_set]))
		down_score = - math.sqrt(len(gene_set_dict[gene_set])/(len(ranked_genes) - len(gene_set_dict[gene_set])))
		running_score = 0
		supremum = 0
		
		for gene in ranked_genes:
			if gene[0] in gene_set_dict[gene_set]:
				running_score += up_score
			else:
				running_score += down_score
			
			if running_score > supremum:
				supremum = running_score
				
		ES.append((gene_set,supremum))
	
	# Sorts the ES list by enrichment score (largest to smallest)
	ES.sort(key=lambda tup: tup[1], reverse=True)

	return ES

def load_KEGG(kegg_file, expression_df):
	"""
		Loads the data from the KEGG gene set file into a dictionary that maps gene set names to genes in that set

		Inputs:
			kegg_file = specially formatted KEGG gene set file
		    expression_df = pandas data frame containing gene expression for each sample
		Returns:
		    gene_set_dict = a dictionary mapping gene_set names to a set of genes in that gene_set
	"""
	gene_set_dict = {}

	with open(kegg_file, 'r') as kegg_f:
		for line in kegg_f:
			line = line.rstrip()
			set_info = line.split('\t')
			genes_in_expression_df = set()
			for gene in set_info[2:]:
				if gene in expression_df['SYMBOL'].unique():
					genes_in_expression_df.add(gene)
			gene_set_dict[set_info[0]] = genes_in_expression_df

	return gene_set_dict

def permutations(expression_df,samples_df,gene_set_dict,num_permute):
	"""
		Conducts num_permute permutations and caculates the enrichment scores for each permutation. The enrichment
		scores from all the permutations are returned in a dictionary

		Inputs:
		   expression_df = pandas data frame containing gene expression for each sample
		   samples_df = pandas data frame containing samples and their classification (1 = patient, 0 = healthy)
		   gene_set_dict = a dictionary mapping gene_set names to a set of genes in that gene_set
		   num_permute = the number of permutations that should be conducted
		Returns:
		   ES_dict = a dictionary mapping gene set names to a list of enrichment scores (one ES for each permutation)
	"""
	permute_samples = samples_df
	ES_dict = {}
	for perumutation in range(num_permute):
		permute_samples[1] = np.random.permutation(permute_samples[1])

		ES = GSEA(expression_df,samples_df,gene_set_dict)

		for set_score in ES:
			if set_score[0] in ES_dict:
				ES_dict[set_score[0]].append(set_score[1])
			else:
				ES_dict[set_score[0]] = [set_score[1]]

	return ES_dict

def print_significant_pathways(ES, ES_dict):
	"""
		Calculates the p-value for each gene set and prints out the number of gene sets that are significant,
		p-value less than 0.05

		Inputs:
		   ES = a list of tuples tha contain the gene set name and the enrichment score of that gene set
		   ES_dict = a dictionary mapping gene set names to a list of enrichment scores (one for each permutation)
		Output:
		   Prints the number of signficant pathways (pathways that have a p-value less than 0.05)
	"""
	p_dict = {}
	for org_score in ES:
		significant_value = 0
		for perm_score in ES_dict[org_score[0]]:
			if org_score[1] <= perm_score:
				significant_value += 1
		p_dict[org_score[0]] = significant_value/100 * len(ES)

	num_pathways = 0
	for p in p_dict:
		if p_dict[p] < 0.05:
			num_pathways += 1
	
	print("Number significant pathways =", num_pathways)

def main():
	# Check that the file is being properly used
	if (len(sys.argv) != 4):
		print("Please specify an expression file, samples file, k, and p as args.")
		return
		
	# Input files
	expression_file = sys.argv[1]
	samples_file = sys.argv[2]
	kegg_file = sys.argv[3]

	# Load files into panda data frames
	expression_df = pd.read_csv(expression_file, sep='\t')
	samples_df = pd.read_csv(samples_file, sep='\t', header=None)

	# Load KEGG file into a dictionary mapping the gene set names to a list of genes in that set
	gene_set_dict = load_KEGG(kegg_file, expression_df)
	
	# Calls the GSEA function to get the original enrichment score
	ES_org = GSEA(expression_df,samples_df,gene_set_dict)

	# Exports the gene set name and the associated enrichment score to a txt file 
	with open('kegg_enrichment_scores.txt', 'w') as output_f:
		for set_score in ES_org:
			output_f.write(set_score[0] + '\t' + str("%.2f" % round(set_score[1],2)) + '\n')

	# Runs 100 permutations of the sample data
	num_permute = 100
	ES_dict = permutations(expression_df,samples_df,gene_set_dict,num_permute)

	# Calculates the p-values and prints the significant pathways 
	print_significant_pathways(ES_org, ES_dict)

if __name__=="__main__":
	main()
