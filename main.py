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

def main():

    targets = ["1b8d_CF","1azs_FD","1d4v_BF"]
    protInfo = textParser()

    train,test = trainTestSplit(protInfo)

    master = {}

    for prot in targets:
        if prot[:-3] not in master:
            master[prot[:-3]] = {}
        chains = prot[-2:]
        if chains[0] not in master[prot[:-3]]:
            master[prot[:-3]][chains[0]] = featureExtract(prot,chains[0])
        if chains[1] not in master[prot[:-3]]:
            master[prot[:-3]][chains[1]] = featureExtract(prot,chains[1])
        print(prot)

    normCols = ["hsed","hseu","seqId","aligns"]
    normaliser = dicNormaliser(columns=normCols,coords=True)
    normaliser.fit(master)

    pd.set_option('display.max_columns', 60)
    master = normaliser.transform(master)

    df = pd.DataFrame()
    for prot in master:
        for chain in master[prot]:
            df = pd.concat([df,master[prot][chain]],axis=0,ignore_index=True)
    print(df.describe())
    return master

    '''
    prefix = "./PPI4DOCK/PPI4DOCK_docking_set/"
    chains = target[-2:]

    chainFile = f"{prefix}/{target}/{chains[0]}_model_st.pdb"

    targetFile = f"{prefix}/{target}/{target}_st.pdb"
    
    msaFile = f"./PPI4DOCK/PPI4DOCK_MSA/{target}/{target[:-2]}{chains[0]}_coMSA.fasta"
    pssm = msa2pssm(msaFile)
    pd.set_option('display.max_columns',61)
    output = rawChainParser(chainFile,chains[0],pssm)

    '''
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