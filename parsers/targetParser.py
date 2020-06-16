from Bio import PDB
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import pi
from .helpers import getSeqIndex, centralCarbon

'''
Functions to parse single target complex into a contact or distance matrix
'''

def distanceParser(targetFile,strippedRes):
    parser = PDB.PDBParser()
    structure = parser.get_structure(targetFile[-14:-3], targetFile)

    model = structure[0]  # only one model
    chainIDs = [chain.id for chain in model]

    dimA = len(model[chainIDs[0]]) + strippedRes[0]
    dimB = len(model[chainIDs[1]]) + strippedRes[1]
    dimensions = (dimA,dimB)

    mat = np.zeros(shape=dimensions)

    for resA in model[chainIDs[0]]:
        for resB in model[chainIDs[1]]:
            mat[getSeqIndex(resA),getSeqIndex(resB)] = residueDistance(resA,resB)

    return mat

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

