"""
This file conducts a molecular dynamics simulation for an input rvc file and outputs an rvc file and an erg file 
every 10 time steps. The output rvc file can be used to visualize the molecular dynamics simulation once it is converted 
to a crd format. The output erg file tracks the bond energy, nonbond energy, kinetic energy, and total energy for every 10 timesteps.

Usage: python MyMD.py --iF input.rvc <other params>

"""

# import numpy as np 
# import argparse

def calcTanimoto(fpA,fpB):
	Tc = size(fpA & fpB)/size(fpA | fpB)
	return round(Tc,6)

def readDrugFile(filename):
	drug_dict = {}
	with open(filename, 'r') as input_f:
		input_f.readline()
		for line in input_f:
			drug = line.strip().split(',')
			drug_dict[drug[0]] = set(drug[2].split())
	
	return drug_dict

def readTargetFile(filename):
	target_dict = {}
	with open(filename, 'r') as input_f:
		input_f.readline()
		for line in input_f:
			target = line.strip().split(',')
			if target[0] in target_dict:
				target_dict[target[0]].add((target[1],target[2]))
			else:
				target_dict[target[0]] = set([(target[1],target[2])])

	return target_dict
