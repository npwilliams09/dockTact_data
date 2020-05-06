from Bio import PDB
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

'''
Functions to parse single target complex into a contact or distance matrix
'''
''' REDUNDANT FUNCTION
def contactParser(targetFile,threshold):
    parser = PDB.PDBParser()
    structure = parser.get_structure(targetFile[-14:-3],targetFile)

    model = structure[0] #only one model
    chainIDs = [chain.id for chain in model]

    dimensions = (len(model[chainIDs[0]]),len(model[chainIDs[1]]))
    mat = np.zeros(shape=dimensions,dtype=bool)


    atoms = PDB.Selection.unfold_entities(model,'A')
    ns = PDB.NeighborSearch(atoms)

    contacts = ns.search_all(threshold,'R')

    i = 0

    for resPair in contacts:
        if resPair[0].get_parent() != resPair[1].get_parent():
            i += 1
            #check chainIDS are in order and reorder if needed
            if resPair[0].get_parent().id != chainIDs[0]:
                resPair = (resPair[1],resPair[0])

            mat[getSeqIndex(resPair[0]),getSeqIndex(resPair[1])] = 1 #set to true for a contact

    f, ax = plt.subplots(figsize=(11, 9))

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(mat)
    plt.show()

    return mat
'''
def distanceParser(targetFile):
    parser = PDB.PDBParser()
    structure = parser.get_structure(targetFile[-14:-3], targetFile)

    model = structure[0]  # only one model
    chainIDs = [chain.id for chain in model]

    dimensions = (len(model[chainIDs[0]]), len(model[chainIDs[1]]))
    mat = np.zeros(shape=dimensions)

    resA = model[chainIDs[0]][1]
    resB = model[chainIDs[1]][2]

    for resA in model[chainIDs[0]]:
        for resB in model[chainIDs[1]]:
            mat[getSeqIndex(resA),getSeqIndex(resB)] = residueDistance(resA,resB)
    sns.heatmap(mat)
    plt.show()

    return mat

def residueDistance(resA,resB):
    return centralCarbon(resA) - centralCarbon(resB)

def centralCarbon(residue):
    if "CA" in residue:
        return residue["CA"]
    else:
        return residue["C"] #special case for glycine residues

def getSeqIndex(residue):
    ids = residue.get_id()
    return ids[1] - 1