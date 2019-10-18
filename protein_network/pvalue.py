"""
This file generates and prints a bootstrapped p-value for a pair of input protiens. The file takes the number of 
bootstrap iterations and the set seed as optional inputs. 

Usage: pvalue.py -n <INT> -r 214 <drugs.csv> <targets.csv> <proteinA> <proteinB>

The  pValue class can also be imported to other files and intialized as follows:

from pvalue import pValue
from argparse import Namespace

pValue(Namespace(
        drugs=<drugs.csv>,
        targets=<targets.csv>,
        proteinA=<proteinA>,
        proteinB=<proteinB>,
        r=214,
        n=<INT>))

"""
import argparse
import numpy as np
from chemoUtils import *

class pValue(object):
    """
    p value class that contains all functions necessary for calculating and 
    outputting a bootstrapped p-value.

    """
    def __init__(self, args):
        """
        Initializes a pValue object, and reads in the drugs.csv and targets.csv files
        using functions in the chemoUtils file.

        Parameters: 
        -----------
        args : a namespace of all the arguments passed in from the command line or another file
        """
        self.args = args
        self.drug_dict = readDrugFile(self.args.drugs)
        self.target_dict = readTargetFile(self.args.targets,'protein:drug')

    def bootstrap(self):
        """
        Conducts a bootstrap for the number of iterations input by the user. The function 
        counts the number of tanimoto bootstrap values are greater than the original 
        tanimoto summary value.

        Returns: 
        --------
        The bootstrapped p-value 
        """
        t_summ = self.calcTSummary(self.target_dict[self.args.proteinA],self.target_dict[self.args.proteinB])
        
        a_size = len(self.target_dict[self.args.proteinA])  # size of ligand set for protein A
        b_size = len(self.target_dict[self.args.proteinB])  # size of ligand set for protein B
        drug_list = list(self.drug_dict.keys())

        boot_count = 0
        np.random.seed(self.args.r) # set random seed
        for i in range(self.args.n):
            a_samp = np.random.choice(drug_list,size=a_size)
            b_samp = np.random.choice(drug_list,size=b_size)

            t_bootstrap = self.calcTSummary(a_samp,b_samp)
            if t_bootstrap >= t_summ:
                boot_count += 1

        return boot_count/self.args.n

    def calcTSummary(self,ligandSetA,ligandSetB):
        """
        Calculates the Tanimoto summary for the input ligand sets. 
        The tanimoto coefficient is calculated using a function from chemoUtils.

        Parameters: 
        -----------
        ligandSetA : set of ligands (size of ligand set for protein A)
        ligandSetB : set of ligands (size of ligand set for protein B)

        Returns: 
        --------
        t_summ : the tanimoto summary value 
        """
        t_summ = 0 
        for a_ligand in ligandSetA:
            for b_ligand in ligandSetB:
                t_c = calcTanimoto(self.drug_dict[a_ligand],self.drug_dict[b_ligand])
                if t_c > 0.5:
                    t_summ += t_c
                
        return t_summ

def parseArgs(**kwargs):
    """
    Uses argparse to read in all the arguments from the command line. For number of iterations and set seed
    variables, the argument is set to a default if not provided. The drugs file, targets file, and two "uniprot_accession" IDs
    are required.
    """
    parser = argparse.ArgumentParser()
    
    # Mandatory inputs
    parser.add_argument("drugs", help="Drugs file", type=str)
    parser.add_argument("targets", help="Targets file", type=str)
    parser.add_argument("proteinA", help="\"uniprot_accession\" ID", type=str)
    parser.add_argument("proteinB", help="\"uniprot_accession\" ID", type=str)

    # Optional variables
    parser.add_argument('-n', type=int, nargs='?', default=500, help='Number of iterations')
    parser.add_argument('-r', type=int, nargs='?', default=214, help='Set seed') 
    
    args = parser.parse_args()
    return pValue(args)

if __name__ == '__main__':
    pFunc = parseArgs()
    print(pFunc.bootstrap())

