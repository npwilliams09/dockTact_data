from Bio import PDB
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import pi

'''
Functions to parse single target complex into a contact or distance matrix
'''

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

def distanceParser(targetFile,strippedRes):
    parser = PDB.PDBParser()
    structure = parser.get_structure(targetFile[-14:-3], targetFile)

    model = structure[0]  # only one model
    chainIDs = [chain.id for chain in model]

    dimA = len(model[chainIDs[0]]) + strippedRes[0]
    dimB = len(model[chainIDs[1]]) + strippedRes[1]
    dimensions = (dimA,dimB)
    # print(dimensions)

    mat = np.zeros(shape=dimensions)

    for resA in model[chainIDs[0]]:
        for resB in model[chainIDs[1]]:
            mat[getSeqIndex(resA),getSeqIndex(resB)] = residueDistance(resA,resB)
    sns.heatmap(mat)


    plt.savefig("distogram.png")

    return mat

def residueDistance(resA,resB):
    return centralCarbon(resA) - centralCarbon(resB)

def centralCarbon(residue):
    if "CA" in residue:
        return residue["CA"]
    else:#special case for glycine residues
        #bipython code sample for virtual beta carbon
        n = residue['N'].get_vector()
        c = residue['C'].get_vector()
        ca = residue['CA'].get_vector()

        #center at origin
        n = n-ca
        c = c-ca
        #find rotation matrix that rotates n
        #-120 degrees along the ca-c vector
        rot = PDB.rotaxis(-pi*120.0/180.0)
        #apply rotation to ca-n vector
        cb_at_origin = n.left_multiply(rot)
        #put on top of ca atom
        cb = cb_at_origin+ca
        return cb


def getSeqIndex(residue):
    ids = residue.get_id()
    return ids[1] - 1

def parseStrippedRes(filePath):
    with open(filePath) as file:
        lineA = stripLine(file.readline())
        lineB = stripLine(file.readline())
        return (lineA,lineB)

def stripLine(line):
    res = []
    for a in [x.strip() for x in line.split("\t")][1:]:
        if (a != "None"):
            ranges = [int(x) for x in a.split('-')]
            res.extend(list(range(ranges[0], ranges[1] + 1)))
    return len(res)

