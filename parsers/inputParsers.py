import numpy as np
from Bio import PDB
from .helpers import getSeqIndex, getResCoords
'''
Functions to parse single inputs
'''

#Returns Lx3 numpy array for input structure
def coordinateParser(inputFile,stripped):
    parser = PDB.PDBParser()
    structure = parser.get_structure(inputFile[-14:-3], inputFile)

    chainID = [chain.id for chain in structure[0]][0]

    chain = structure[0][chainID]

    seqLength = len(chain) + stripped

    output = np.zeros(shape=(seqLength,3))

    for residue in chain:
        coords = getResCoords(residue).reshape(1,3)
        output[getSeqIndex(residue),:] = coords
    output.reshape((1,-1,3))
    return output

