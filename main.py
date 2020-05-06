from parsers.textParser import textParser
from parsers.targetParser import contactParser, distanceParser

def main():
    threshold = 15
    #contactParser("./PPI4DOCK/PPI4DOCK_docking_set/1a2y_AC/1a2y_AC_st.pdb",threshold)
    distanceParser("./PPI4DOCK/PPI4DOCK_docking_set/1a2y_AC/1a2y_AC_st.pdb")


if __name__ == '__main__':
    main()