"""
Protein Sequencing Project
Name:
Roll Number:
"""

from dataclasses import replace
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file = open(filename,"r").read()
    stri = ""
    for line in file.splitlines():
        stri = stri+line
    return stri

'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    condonslst = []
    var = ["UGA","UAG","UAA"]
    for word in range(startIndex,len(dna),3):
        dna = dna.replace("T","U")
        condns = dna[word:word+3]
        if condns not in var:
            condonslst.append(condns)
        else:
            condonslst.append(condns)
            break
    return condonslst
'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    file = open(filename,"r")
    protns = json.load(file)
    codondictnry= {}
    for key in protns:
        for values in protns[key]:
            values = values.replace("T","U")
            codondictnry[values] = key
    return codondictnry

'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''

def generateProtein(codons, codonD):
    protnlst = []
    for rna in codons:
        for rnaproteins in codonD:
            if rna == rnaproteins:
                protnlst.append(codonD[rnaproteins])
                if protnlst[0] == "Met":
                    protnlst[0] = "Start"
    return protnlst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    rnalst = readFile(dnaFilename)
    protenlst = makeCodonDictionary(codonFilename)
    totalprotenlst = []
    i = 0
    unusedltrs = 0
    while i<len(rnalst) :
        word = rnalst[i:i+3]
        if word == "ATG":
            rna = dnaToRna(rnalst, i)
            totalprotenlst.append(generateProtein(rna,protenlst))
            i = i+3*len(rna)
        else:
            unusedltrs += 1
            i += 1      
    return totalprotenlst


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    cmnprotenlst = []
    for proten1 in proteinList1:
        for proten2 in proteinList2:
            if proten1 == proten2 and proten1 not in cmnprotenlst:
                cmnprotenlst.append(proten1)
    return cmnprotenlst


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    protenlst = []
    for lst in proteinList:
        for protens in lst:
            lst.append(protens)
    return protenlst


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoaciddictnry = {}
    for word in aaList:
        if word not in aminoaciddictnry:
            aminoaciddictnry[word] = 0
        if word in aminoaciddictnry:
            aminoaciddictnry[word] += 1
    return aminoaciddictnry


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    return


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()

    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """
    test.testCommonProteins()
    # test.testCombineProteins()
    test.testAminoAcidDictionary()


    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
