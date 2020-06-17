import numpy as np
from Bio import PDB,AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.PDB.ResidueDepth import ResidueDepth
import pandas as pd
from .helpers import getSeqIndex, getResCoords, centralCarbon
from sklearn.preprocessing import LabelBinarizer
'''
Functions to parse single inputs
'''

#Returns a dataframe with raw data for a single chain
def rawChainParser(filepath, chainID, pssm):
    parser = PDB.PDBParser()
    structure = parser.get_structure(chainID, filepath)
    model = structure[0]

    d = PDB.DSSP(model, filepath)

    rd = ResidueDepth(model)

    hs = PDB.HSExposure.HSExposureCA(model)

    df = pd.DataFrame()
    acids = getAminoAcids(structure)

    for i,residue in enumerate(model[chainID]):
        seqID = getSeqIndex(residue)

        row = {}
        x,y,z = getResCoords(residue)
        row["AA"] = acids[i]
        row["x"] = x[0]
        row["y"] = y[0]
        row["z"] = z[0]

        tupKey = (chainID, (' ', seqID, ' '))
        row["res_depth"] = rd[tupKey][0]
        row["ca_depth"] = rd[tupKey][1]

        dssp = d[(chainID,seqID)]

        row["ss"] = dssp[2]
        row["asa"] = dssp[3]
        row["phi"] = dssp[4]/360.0
        row["psi"] = dssp[5]/360.0

        if tupKey in hs:
            row["hseu"] = hs[tupKey][0]
            row["hsed"] = hs[tupKey][1]
        else: #NO HSE CALCULATED
            row["hseu"] = 0
            row["hsed"] = 0

        row["seqId"] = seqID -1
        #row["bFactor"] = centralAtom.get_bfactor() #must be done using zhang lab tool resQ instead

        #get pssm row
        pssmRow = pssm[seqID-1]

        row["aligns"] = sum(pssmRow.values())  # total alignments
        if(row["aligns"]!= 0):
            for key in pssmRow.keys():
                if key not in list('ABCDEFGHIKLMNPQRSTVWYZ'):
                    print(key)
                row["pssm_" + key] = pssmRow[key]/row["aligns"]



        if (residue.is_disordered()):
            print(f"disorded atom in res {getSeqIndex()}")

        df = df.append(row,ignore_index=True)

    #check for missed columns in pssm
    for a in 'ABCDEFGHIKLMNPQRSTVWYZ':  # check for missing aas
        if ("pssm_" + a) not in df:
            df["pssm_" + a] = 0.0

    aaEncoder = AAonehot()
    aaTransformed = aaEncoder.transform(df["AA"])
    aaCols = ["AA_"+x for x in aaEncoder.classes_]
    aaDf = pd.DataFrame(aaTransformed,columns=aaCols)

    ssEncoder = ssOneHot()
    ssTransformed = ssEncoder.transform(df["ss"])
    ssCols = ["SS_" + x for x in ssEncoder.classes_]
    ssDf = pd.DataFrame(ssTransformed,columns=ssCols)

    #center coordinates
    max = df.max()
    min = df.min()
    df["x"] = df["x"] - (max["x"] + min["x"]) / 2
    df["y"] = df["y"] - (max["y"] + min["y"]) / 2
    df["z"] = df["z"] - (max["z"] + min["z"]) / 2

    df = pd.concat([df,aaDf,ssDf],axis=1).drop(["AA","ss"],axis=1)

    if(df.shape[1] != 62):
        print(df)
    assert df.shape[1] == 62, f"Incorrect pssmdf shape = {df.shape[1]} for file: {filepath}" #error check
    return df

#returns df of pssm given a fasta MSA as input
def msa2pssm(msaFile):
    alignment = AlignIO.read(msaFile, "fasta")
    summary = AlignInfo.SummaryInfo(alignment)
    return summary.pos_specific_score_matrix(chars_to_ignore=["-","X"])

#returns sequence as a Lx1 numpy array
def getAminoAcids(structure):
    ppb = PDB.PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(structure):
        sequence += str(pp.get_sequence())
    return list(sequence)

#returns the onehot encoder for amino acids
def AAonehot():
    AAstr = 'ACDEFGHIKLMNPQRSTVWY'
    AAs = pd.Series(list(AAstr))
    return LabelBinarizer().fit(AAs)

#returns the onehot encoder for secondary structure
def ssOneHot():
    ssCats = ['H','B','E','G','I','T','S','-']
    ssCats = pd.Series(ssCats)
    return LabelBinarizer().fit(ssCats)
