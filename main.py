from parsers.textParser import textParser
from parsers.targetParser import contactParser, distanceParser, parseStrippedRes

def main():
    target = "1axi_BD"
    prefix = "./PPI4DOCK/PPI4DOCK_docking_set/"
    strippedRes = parseStrippedRes(prefix + f"{target}/stripped_res.txt")
    distanceParser(prefix + f"{target}/{target}_st.pdb",strippedRes)



if __name__ == '__main__':
    main()