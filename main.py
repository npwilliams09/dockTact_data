from parsers.textParser import textParser
from parsers.targetParser import distanceParser, parseStrippedRes
from parsers.inputParsers import rawChainParser, getAminoAcids, msa2pssm
from parsers.helpers import getSeqIndex, getResCoords
from utilities.dataSplit import trainTestSplit
from utilities.dicNormaliser import dicNormaliser
import os
from Bio import PDB, SeqIO,AlignIO

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gc

def main():
    protInfo = textParser()

    train,test = trainTestSplit(protInfo)

    trainDic = loadSet(train,breakPoint=5)

    normCols = ["hsed","hseu","seqId","aligns"]
    normaliser = dicNormaliser(columns=normCols,coords=True)

    trainDic = normaliser.fit_transform(trainDic)
    saveSet(trainDic)
    del trainDic
    gc.collect()

    testDic = loadSet(test,breakPoint=5)
    testDic = normaliser.transform(testDic)
    saveSet(testDic)
    del testDic
    gc.collect()

def loadSet(setList,breakPoint=10000):
    master = {}
    entries = max(len(setList),breakPoint)
    i = 0

    for prot in setList:
        if prot[:-3] not in master:
            master[prot[:-3]] = {}
        chains = prot[-2:]
        if chains[0] not in master[prot[:-3]]:
            master[prot[:-3]][chains[0]] = featureExtract(prot, chains[0])
        if chains[1] not in master[prot[:-3]]:
            master[prot[:-3]][chains[1]] = featureExtract(prot, chains[1])

        i+=1
        if (i>breakPoint):
            return master
        print(f'\r>>Loading Set: {"{:.2f}".format(i/entries)}%')

    return master

def loadTargets(complexes):


def saveSet(dic):
    entries = len(dic)
    i = 0

    makeFolder("./output/chains")
    makeFolder("./output/targets")

    for prot in dic:
        for chain in dic[prot]:
            filename = f"{prot}_{chain}.npy"
            with open("./output/chains/"+filename,'wb') as f:
                np.save(f,dic[prot][chain].values)

def featureExtract(target, chain):
    prefix = "./PPI4DOCK/PPI4DOCK_docking_set/"
    chainFile = f"{prefix}/{target}/{chain}_model_st.pdb"

    msaFile = f"./PPI4DOCK/PPI4DOCK_MSA/{target}/{target[:-2]}{chain}_coMSA.fasta"
    pssm = msa2pssm(msaFile)

    output = rawChainParser(chainFile, chain, pssm)
    return output

def makeFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)

if __name__ == '__main__':
    main()