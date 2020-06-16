import json
'''
parses the "PPI4DOCK_list.txt" file into a .json file for convenience
and returns the output dictionary
'''
def textParser():
    fileAddr = "./PPI4DOCK/PPI4DOCK_list.txt"
    dic = {}

    file = open(fileAddr,"r")
    headers = file.readline().split(" ")[1:] #extract header titles
    headers = [(x.split(".",1)[-1]).strip() for x in headers] #remover numbering and new lines

    for line in file: #line by line now header removed
        data = line.split(" ") #split into list
        for i,variable in enumerate(data):
            if (i==0): #create the sub dictionary
                dic[data[0]] = {}
            dic[data[0]][headers[i]] = variable.strip() #add data next to each relevant header

    with open("./output/proteinInfo.json",'w') as writer: #write output to json file
        json.dump(dic,writer,ensure_ascii=False, indent=4)

    return dic #return the output dictionary
