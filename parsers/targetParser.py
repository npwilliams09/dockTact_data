from Bio import PDB
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import pi
from .helpers import getSeqIndex, centralCarbon

'''
Functions to parse single target complex into a contact or distance matrix
'''

def distanceParser(targetFile,chains, seqIDs):
    parser = PDB.PDBParser()
    structure = parser.get_structure(targetFile[-14:-3], targetFile)

    model = structure[0]  # only one model

    dimA = len(seqIDs[0])
    dimB = len(seqIDs[1])
    dimensions = (dimA,dimB)

    mat = np.zeros(shape=dimensions)

    for i,idA in enumerate(seqIDs[0]):
        for j,idB in enumerate(seqIDs[1]):
            resA = model[chains[0]][idA]
            resB = model[chains[1]][idB]
            mat[i,j] = residueDistance(resA,resB)
    return mat

def dis2contact(mat,threshold=8):
    return np.where(mat <= threshold, 1, 0)

def residueDistance(resA,resB):
    return centralCarbon(resA) - centralCarbon(resB)

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

