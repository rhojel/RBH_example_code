"""
Generates tanimoto coefficients for ever pair of drugs in the drugs.csv. 
Outputs the drug pairsm the tanimoto coefficient, and whether or not the drugs have a 
target in common into an output file.

Ex. Output File
----------------
DB00001,DB00006,0.810127,1 
DB00001,DB00007,0.835443,0 
DB00001,DB00010,0.691358,0

Usage: tanimoto.py <drugs.csv> <targets.csv> <outputfile.csv>

"""
from chemoUtils import *
import sys

class tanimotoGeneration(object):
    """
    Tanimoto generation class that contains all functions necessary for calculating and 
    outputting the tanimoto file.

    """
    def __init__(self, drugs, targets):
        """
        Initializes a tanimotoGeneration object and reads in the drugs.csv and targets.csv files
        using functions in the chemoUtils file.

        Parameters: 
        -----------
        drugs : the directory of the input file (drugs.csv)
        targets : the directory of the input file (targets.csv)
        """
        self.drug_dict = readDrugFile(drugs)
        self.target_dict = readTargetFile(targets,'drug:protein')

    def fileOuput(self,output_file):
        """
        Calculates the tanimoto ceofficients for each pair and outputs the tanimoto file.
        The tanimoto coefficient is calculated using a function from chemoUtils.

        Parameters: 
        -----------
        output_file : the name and directory of the output file

        """
        completedPairs = set()  # Tracks previously seen drug pairs

        with open(output_file, 'w') as output_f:
            for drug in self.drug_dict:
                for pair in self.drug_dict:
                    sharedTarget = 0
                    if (pair,drug) not in completedPairs and drug != pair:
                        Tc = calcTanimoto(self.drug_dict[drug],self.drug_dict[pair])
                        if drug in self.target_dict and pair in self.target_dict:
                            if len(self.target_dict[drug] & self.target_dict[pair]) > 0:
                                sharedTarget = 1
                        output_f.write(str(drug) +','+str(pair)+','+'%.6f'%Tc+','+str(sharedTarget)+'\n')
                        completedPairs.add((drug,pair))

def main():
    # Check that the file is being properly used
    if (len(sys.argv) != 4):
        print("Please specify drugs.csv, targets.csv, and output file name.")
        return
    
    tG = tanimotoGeneration(sys.argv[1],sys.argv[2])
    tG.fileOuput(sys.argv[3])

if __name__ == '__main__':
    main()