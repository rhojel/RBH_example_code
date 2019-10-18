#	Implementation of the A-Priori frequent item set algorithm to identify frequent pairs and frequent triples 
#	Exports files with association rules of frequent pairs/frequent triples and their confidence values

#	By Richard Hojel
#	For CS 246 Homework 1 Problem 2
import sys

def calcItemFrequency(input_file, s):
	freqItems = {}
	with open(input_file, 'r') as input_f:
		for line in input_f:
			basket = line.split()
			for item in basket:
				if item in freqItems:
					freqItems[item] += 1
				else:
					freqItems[item] = 1

	for item, value in list(freqItems.items()):
		if value < s:
			del freqItems[item]

	return freqItems

def freqPairSet(input_file, s, freqItems):
	freqPairs = {}
	with open(input_file, 'r') as input_f:
		for line in input_f:
			basket = line.split()
			for item1 in basket:
				if item1 in freqItems:
					for item2 in basket:
						if item2 in freqItems and item1 < item2:
							if (item1,item2) in freqPairs:
								freqPairs[(item1,item2)] += 1
							else:
								freqPairs[(item1,item2)] = 1

	for pair, value in list(freqPairs.items()):
		if value < s:
			del freqPairs[pair]

	return freqPairs

def freqTripleSet(input_file, s, freqPairs):
	freqTriples = {}
	with open(input_file, 'r') as input_f:
		for line in input_f:
			basket = line.split()
			for item1 in basket:
					for item2 in basket:
						if (item1,item2) in freqPairs:
							for item3 in basket:
								if (item1,item3) in freqPairs and (item2,item3) in freqPairs and (item1<item2<item3):
									if (item1,item2,item3) in freqTriples:
										freqTriples[(item1,item2,item3)] += 1
									else:
										freqTriples[(item1,item2,item3)] = 1

	for triple, value in list(freqTriples.items()):
		if value < s:
			del freqTriples[triple]

	return freqTriples

def pairConfidence(freqPairs, freqItems):
	# Calculates confidence of all possible association rules of pairs in dict 'freqPairs'
	confPairs = {}
	for pair in freqPairs:
		confPairs[pair] = freqPairs[pair]/freqItems[pair[0]]
		confPairs[(pair[1],pair[0])] = freqPairs[pair]/freqItems[pair[1]]
	return confPairs

def tripleConfidence(freqTriples, freqPairs):
	# Calculates confidence of all possible association rules of triples in dict 'freqTriples'
	confTriples = {}
	for triple in freqTriples:
		confTriples[triple] = freqTriples[triple]/freqPairs[(triple[0],triple[1])]
		confTriples[(triple[0],triple[2],triple[1])] = freqTriples[triple]/freqPairs[(triple[0],triple[2])]
		confTriples[(triple[1],triple[2],triple[0])] = freqTriples[triple]/freqPairs[(triple[1],triple[2])]
	return confTriples

def main():
	# Read in file
	input_file = sys.argv[1]
	s = 100

	# A-priori first pass, calculate the frequent items (frequency above threshold 's')
	freqItems = calcItemFrequency(input_file, s)

	# A-priori second pass, calculate the frequent item pairs (frequency above threshold 's')
	# Only considers candidate pairs that are both in the frequent item set
	freqPairs = freqPairSet(input_file, s, freqItems)

	# A-priori third pass, calculate the frequent item triples (frequency above threshold 's')
	# Only considers candidate triples for which all the pairs are in the frequent pairs set
	freqTriples = freqTripleSet(input_file, s, freqPairs)
	
	# Problem 2(d)
	# Calculates the confidence of all possible association rules in the frequent pairs set
	# Dict format- (a,b):calcConfidence where calcConfidence is the confidence for the association rule: a -> b
	confPairs = pairConfidence(freqPairs, freqItems)

	# Problem 2(e)
	# Calculates the confidence of all possible association rules in the frequent triples set
	# Dict format- (a,b,c):calcConfidence where calcConfidence is the confidence for the association rule: (a,b) -> c
	confTriples = tripleConfidence(freqTriples, freqPairs)

	# UNCOMMENT to write output Files
	# Writes the two dictionaries of association rules and confidence values to a file
	with open('confPairs.txt', 'w') as out_file:
		for key, value in confPairs.items():
			out_file.write(str(key[0])+'->'+str(key[1])+'\t'+str(value)+'\n')

	with open('confTriples.txt', 'w') as out_file:
		for key, value in confTriples.items():
			out_file.write('('+str(key[0])+','+str(key[1])+')'+'->'+str(key[2])+'\t'+str(value)+'\n')		

if __name__ == '__main__':
	main()