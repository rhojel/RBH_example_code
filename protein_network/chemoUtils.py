"""
This file contains helper functions used by tanimoto.py, pvalue.py, networkgen.py, and plot_graph.py.

The first function returns the tanimoto coefficient and the remaining functions import a file 
and return a data structure containing important file information.

"""
def calcTanimoto(fpA,fpB):
    """
    Computes the tanimoto coefficient (Jaccard Index) for the fingerprints of two molecules.

    Parameters: 
    -----------
    fpA : fingerprint (set of keys indicating the presence of a feature) of molecule A 
    fpB : fingerprint (set of keys indicating the presence of a feature) of molecule B

    Returns: 
    --------
    The tanimoto coefficient rounded to six decimal places.
    """
    Tc = len(fpA & fpB)/len(fpA | fpB)
    return round(Tc,6)

def readDrugFile(filename):
    """
    Reads the drug.csv file into a dictionary.

    Parameters: 
    -----------
    filename : the directory of the input file (drugs.csv)

    Returns: 
    --------
    drug_dict : a dictionary mapping a drug to its fingerprint (drug:fingerprint)
    """
    drug_dict = {}
    with open(filename, 'r') as input_f:
        input_f.readline()
        for line in input_f:
            drug = line.strip().split(',')
            drug_dict[drug[0]] = set(drug[2].split())
    
    return drug_dict

def readTargetFile(filename, form):
    """
    Reads the targets.csv file into a dictionary.

    Parameters: 
    -----------
    filename : the directory of the input file (targets.csv)
    form : indicates the organization of the dictionary (either drug:protein or protein:drug)

    Returns: 
    --------
    target_dict : a dictionary mapping either a drug to a set of proteins (uniprot_accession) 
    or a protien (uniprot_accession) to a set of drugs
    """
    index = 1
    if form == 'drug:protein':
        index = 0

    target_dict = {}
    with open(filename, 'r') as input_f:
        input_f.readline()
        for line in input_f:
            target = line.strip().split(',')
            if target[index] in target_dict:
                target_dict[target[index]].add(target[abs(index-1)])
            else:
                target_dict[target[index]] = set([target[abs(index-1)]])

    return target_dict

def readNodeFile(filename, colorIndicator = False):
    """
    Reads the protein_nodes.csv file into one or two dictionaries. 

    Parameters: 
    -----------
    filename : the directory of the input file (protein_nodes.csv)
    colorIndicator : indicates whether or not to make and return a color dictionary that maps uniprot_id:indications

    Returns: 
    --------
    nodes : a dictionary mapping uniprot_accession to uniprot_id
    color : a dictionary mapping uniprot_id to indications
    """
    color = {}
    nodes = {}
    with open(filename, 'r') as input_f:
        input_f.readline()
        for line in input_f:
            node = line.strip().split(',')
            nodes[node[0]] = node[1]
            if colorIndicator: 
                color[node[1]] = node[2]
    if colorIndicator:
        return nodes, color 
    else:
        return nodes

