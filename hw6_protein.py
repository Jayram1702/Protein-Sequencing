"""
Protein Sequencing Project
Name:
Roll Number:
"""

from dataclasses import replace
from os import remove

from matplotlib.cbook import file_requires_unicode
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
            protenlst.append(protens)
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
    amnofreqlst =[]
    totlcount1 = len(combineProteins(proteinList1))
    totlcount2 = len(combineProteins(proteinList2))
    word = combineProteins(proteinList1) 
    word1 = combineProteins(proteinList2)
    count = aminoAcidDictionary(word)
    count1 = aminoAcidDictionary(word1)
    totalwrd = word + word1
    totalwrd = list(set(totalwrd))
    for i in totalwrd:
        freq1 = 0 
        freq2 = 0
        if i != "Start" and i != "Stop":
            if i not in amnofreqlst:
                if i in word and i in word1:
                    freq1 = count[i]/totlcount1
                    freq2 = count1[i]/totlcount2
                elif i in word and i not in word1:
                    freq1 = count[i]/totlcount1
                    freq2 = 0
                elif i not in word and i in word1:
                    freq1 = 0 
                    freq2 = count1[i]/totlcount2
                freqdiff = freq1-freq2
                if freqdiff > cutoff or freqdiff < -cutoff:
                    amnofreqlst.append([i,freq1,freq2])
    return amnofreqlst


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    lst = []
    lst1 = []
    strng = ''
    print("The following proteins occurred in both DNA Sequences:")
    for i in commonalities:
        if i not in lst:
            lst.append(i[1:-1])
    for j in lst:
        if len(j) > 0:
            j = '-'.join(j)
            lst1.append(j)
            lst1.sort()
    for k in lst1:
        strng += ' '+ k +"\n"
    print(strng)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for a in differences:
        wrd = a[0]
        seq1 = round(a[1]*100,2)
        seq2 = round(a[2]*100,2)
        print(f"{wrd}:{seq1}% in Seq1, {seq2}% in seq2")
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
    labels = []
    wrd = combineProteins(proteinList1) + combineProteins(proteinList2)
    for i in wrd:
        if i not in labels:
            labels.append(i)
    totallabels = sorted(labels)
    return totallabels


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    freqlst = []
    wrd = combineProteins(proteinList)
    totlcnt = len(combineProteins(proteinList))
    cnt = aminoAcidDictionary(wrd)
    for i in labels:
        if i in cnt:
            freq = cnt[i]/totlcnt
            freqlst.append(freq)
        else:
            freqlst.append(0)
    return freqlst


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w = 0.35
    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1, edgecolor = edgeList)
    plt.bar(xLabels, freqList2, width= w, align='edge', label=label2, edgecolor = edgeList)

    plt.xticks(rotation="vertical")
    plt.legend()
    plt.title("Compare Human and Elephant Genes")

    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    colorlst = []
    wordlst1 = []
    for word1 in labels:
        for word in biggestDiffs:
            wordlst1.append(word[0])
        if word1 in wordlst1:
            colorlst.append("black")
        else:
            colorlst.append("white")
    return colorlst

'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    protns1 = synthesizeProteins("data/Human_p53.txt", "data/codon_table.json")
    protns2 = synthesizeProteins("data/Elephant_p53.txt", "data/codon_table.json")
    commonprotns = commonProteins(protns1, protns2)
    amodiff = findAminoAcidDifferences(protns1,protns2,0.005)
    displytext = displayTextResults(commonprotns, amodiff)
    amolabels = makeAminoAcidLabels(protns1,protns2)
    f1 = setupChartData(amolabels, protns1)
    f2 = setupChartData(amolabels, protns2)
    edges = makeEdgeList(amolabels, amodiff)
    finalchart = createChart(amolabels, f1, 'Human', f2, 'Elephant', edgeList=edges)
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
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()

    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()


    ## Uncomment these for Week 3 ##

    # print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    # test.week3Tests()
    # print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
    # test. testMakeAminoAcidLabels()
    # test.testSetupChartData()
    # test.testCreateChart()
    # test.testMakeEdgeList()