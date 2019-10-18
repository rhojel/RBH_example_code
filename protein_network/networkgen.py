"""
This file generates an edge list of significant connections between two proteins (from the protein__nodes.csv file). 
The significant connections are determined using the bootstrapped p-value generated in the pvalue.py file.

Usage: networkgen.py <drugs.csv> <targets.csv> <protein_nodes.csv> 

"""
from pvalue import pValue
from chemoUtils import readNodeFile
from argparse import Namespace
import sys

def createEdgeList(drugs,targets,nodes):
    """
    Creates a set of tuples with each tuple containing two proteins that have an edge between them.
    The protein pairs are deterimined using a bootstrapped pValues (only includes p-values < 0.05). 

    Parameters: 
    -----------
    drugs : the directory of the input file (drugs.csv)
    targets : the directory of the input file (targets.csv)
    nodes : a dictionary that maps uniprot_accession to uniprot_id (only keys are used)

    Returns: 
    --------
    edgeList : a set of tuples, each tuple represents two proteins that have an edge between them 
    """
    seenNodes = set()
    edgeList = set()
    for node in nodes:
        for pair in nodes:
            if node != pair and (pair,node) not in seenNodes:
                pCalcFunc = pValue(Namespace(drugs=drugs,targets=targets,proteinA=node,proteinB=pair,r=214,n=500))
                p = pCalcFunc.bootstrap()
                if p <= 0.05:
                    edgeList.update([(node,pair)])
                
                seenNodes.update([(node,pair)])
    return edgeList

def exportEdgeList(edgeList):
    """
    Outputs the edge list to a txt file

    Parameters: 
    -----------
    edgeList : a set of tuples, each tuple represents two proteins that have an edge between them 
    """
    with open('network_edgelist.txt', 'w') as output_f:
        for edge in edgeList:
            output_f.write(edge[0]+' '+edge[1]+'\n')

def main():
    # Check that the file is being properly used
    if (len(sys.argv) != 4):
        print("Please specify drugs.csv, targets.csv, and protein_nodes.csv.")
        return

    nodes = readNodeFile(sys.argv[3])   # function from chemoUtils file
    edgeList = createEdgeList(sys.argv[1],sys.argv[2],nodes)
    exportEdgeList(edgeList)

if __name__ == '__main__':
    main()
