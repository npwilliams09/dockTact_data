#Retrieves the sequence number of a residue
def getSeqIndex(residue):
    ids = residue.get_id()
    return ids[1]

# Returns residue coords as a tuple in x,y,z format
def getResCoords(res):
    atom = centralCarbon(res)
    coords = atom.get_coord()
    coords = coords.reshape(1,3).T
    return coords[0],coords[1],coords[2]

#Returns the central atom of a residue
def centralCarbon(residue):
    if "CB" in residue:
        return residue["CB"]
    elif "CA" in residue:#special case for glycine residues
        #return alpha carbon
        return residue["CA"]
    else:
        for atom in residue:
            print(f"Error in residue {getSeqIndex(residue)}: {residue}")
            return atom # return first atom for corrupted files

def getPharmacophoreDict():
    fingerprint = {}
    fingerprint["A"] = [1, 0, 0, 1, 1, 0, 0, 2, 13]
    fingerprint["V"] = [3, 0, 0, 1, 1, 0, 0, 2, 19]
    fingerprint["L"] = [4, 0, 0, 1, 1, 0, 0, 2, 22]
    fingerprint["G"] = [0, 0, 0, 1, 1, 0, 0, 2, 10]
    fingerprint["S"] = [0, 0, 0, 1, 2, 0, 0, 3, 14]
    fingerprint["W"] = [10, 0, 0, 1, 2, 9, 0, 2, 27]
    fingerprint["T"] = [1, 0, 0, 1, 2, 0, 0, 3, 17]
    fingerprint["Q"] = [2, 0, 0, 2, 2, 0, 0, 3, 20]
    fingerprint["E"] = [2, 0, 2, 3, 1, 0, 0, 3, 19]
    fingerprint["C"] = [1, 0, 0, 1, 1, 0, 1, 3, 14]
    fingerprint["R"] = [3, 2, 0, 1, 4, 0, 0, 3, 26]
    fingerprint["P"] = [3, 0, 0, 1, 1, 0, 0, 2, 17]
    fingerprint["D"] = [1, 0, 2, 3, 1, 0, 0, 3, 16]
    fingerprint["F"] = [7, 0, 0, 1, 1, 6, 0, 2, 23]
    fingerprint["I"] = [4, 0, 0, 1, 1, 0, 0, 2, 22]
    fingerprint["H"] = [4, 2, 0, 1, 3, 5, 0, 2, 20]
    fingerprint["N"] = [1, 0, 0, 2, 3, 0, 0, 3, 17]
    fingerprint["M"] = [3, 0, 0, 1, 1, 0, 1, 2, 20]
    fingerprint["Y"] = [7, 0, 0, 1, 2, 6, 0, 2, 24]
    fingerprint["K"] = [3, 1, 0, 1, 2, 0, 0, 3, 24]
    # Hydrophobics,Positives,Negatives,Acceptors,Donnors,Aromatics,Sulfurs,Neutral,No. of atoms

    acidDic = {}
    for a in fingerprint.keys():
        dic = {}
        array = fingerprint[a]
        total = array[8]
        dic["Hydrophobic"] = array[0]/total
        dic["Positives"] = array[1]/total
        dic["Negatives"] = array[2]/total
        dic["Acceptors"] = array[3]/total
        dic["Donors"] = array[4]/total
        dic["Aromatics"] = array[5]/total
        dic["Sulfurs"] = array[6]/total
        dic["Neutral"] = array[7]/total
        dic["Total_Atoms"] = total/27 #27 is highest number of atoms, therefore normalises data
        acidDic[a] = dic
    return acidDic