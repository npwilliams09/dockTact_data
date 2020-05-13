def getSeqIndex(residue):
    ids = residue.get_id()
    return ids[1] - 1

def getResCoords(res):
    atom = centralCarbon(res)
    coords = atom.get_coord()
    coords = coords.reshape(1,3).T
    return coords

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