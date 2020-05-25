
# method to be used
def CutLastField(istr):
    cutto = istr.rfind('_')
    return istr[0:cutto]


def GrepStrainFromAlignment(alnfile, strain, afixCut) :
    result = ["", 0, True]

    isTargetTitle = {}
    isTargetTitle[strain] = True

    fastaReader = open(alnfile, "r")
    seq = []
    line = fastaReader.readline()
    while (line != ''):
        if(line.startswith(">")):
            seqtitle = line.strip()[1:].split(" ")[0]
            if afixCut:
                seqtitle = CutLastField(seqtitle)
            seq = []
            line = fastaReader.readline()
            while (line != ''):
                seq.append(line.strip())
                line = fastaReader.readline()
                if(line.startswith(">")):
                    break
            
            sequence = "".join(seq)
            result[1] = len(sequence)
            if(seqtitle in isTargetTitle):
                result[0] = sequence
                result[2] = False
                break
        else :
            line = fastaReader.readline()
    fastaReader.close()

    if result[2]:
        result[0] = '-' * result[1]

    return result

def LimitNucltodiedCodes(ntseq):
    validBase = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'N', 'n', '-']
    ntlist = list(ntseq)
    for idx in range(len(ntlist)):
        if ntlist[idx] not in validBase:
            ntlist[idx] = 'N'
    return ''.join(ntlist)



# parse input
import argparse

parser = argparse.ArgumentParser(description = "Input specifications ")
parser.add_argument("--seqfiles", dest="seqfileListFile", required=True, help="single column text file containing the list of sequence alignment file paths", type=str)
parser.add_argument("--strains", dest="strainListFile", required=True, help="single column text file containing the list of strains to be included in the output", type=str)
parser.add_argument("--afixcut", dest="needOrfAfixCut", required=False, help="--afixcut if sequence IDs in the alignment files are ORF IDs [with afix _1 _2 _3...] rather than the genome IDs > then we'll cut the afix out", default = False, action = "store_true")
parser.add_argument("--xmfa", dest="outputInXmfa", required=False, help="--xmfa if you want block-by-block XMFA format output instead of concatenated supermatrix", default=False, action="store_true")
parser.add_argument("--noIUPAC", dest="iupacNotAllowed", required=False, help="--noIUPAC if you want only [ATCGatcgNn-] and nothing else in your sequences; otherwise by default we leave the sequence as is", default = False, action = "store_true")
parser.add_argument("--out", dest="outputFile", required=True, help="the output concatenated fasta", type=str)

args = parser.parse_args()

seqfileListFile = args.seqfileListFile
strainListFile = args.strainListFile
outputFile = args.outputFile
needOrfAfixCut = args.needOrfAfixCut
outputInXmfa = args.outputInXmfa
iupacNotAllowed = args.iupacNotAllowed

if outputInXmfa:
    print(" use XMFA output format")

if needOrfAfixCut:
    print(" ORF _index part will be cut out")

if iupacNotAllowed:
    print(" Only allow ATCGatcgNn- in the sequences, so we'll convert anything else to N")


# create a string:integer dictionary and an integer:string (reverse)dictionary of strain list
strainListIndex = int
strainListIndex = 0
dictStrainIndex = {}
dictIndexStrain = {}
listReader = open(strainListFile, "r")
for line in listReader:
    strain = line.strip()
    dictStrainIndex[strain] = strainListIndex
    dictIndexStrain[strainListIndex] = strain
    strainListIndex += 1
print("strain list size = " + str(len(dictIndexStrain)))
listReader.close()


# create a list of seq files
seqfileList = []
listReader = open(seqfileListFile, "r")
for line in listReader:
    seqfile = line.strip()
    seqfileList.append(seqfile)
print("alignment file list size = " + str(len(seqfileList)))
listReader.close()


# do concatenation
if outputInXmfa:
    # XMFA case: simply write the same fasta into a block, except that you put gap string '-' * length for the missing strains
    outputWriter = open(outputFile, "w")
    for seqfile in seqfileList:
        for strainIndex in dictIndexStrain:
            strainName = dictIndexStrain[strainIndex]
            strainSeq = GrepStrainFromAlignment(seqfile, strainName, needOrfAfixCut)
            strainSeqNt = ""
            if iupacNotAllowed: 
                strainSeqNt = LimitNucltodiedCodes(strainSeq[0])
            else:
                strainSeqNt = strainSeq[0]
            outputWriter.write(">" + strainName + "\n" + strainSeqNt + "\n")
        outputWriter.write("= genefamily=" + seqfile + "\n")
    outputWriter.close()

else:
    # Standard case:
    # for each strain, collect the sequence and write the lines to the output
    outputWriter = open(outputFile, "w")
    for strainIndex in dictIndexStrain:
        strainName = dictIndexStrain[strainIndex]
        totalLength = 0
        numMissing = 0

        outputWriter.write(">" + strainName + "\n")
        for seqfile in seqfileList:
            concatee = GrepStrainFromAlignment(seqfile, strainName, needOrfAfixCut) # the returning list from the method is [{string, sequence to concatenate}, {integer, length}, {boolean, missing}]
            concateeNt = ""
            if iupacNotAllowed: 
                concateeNt = LimitNucltodiedCodes(concatee[0])
            else:
                concateeNt = concatee[0]
            outputWriter.write(concateeNt)
            totalLength += concatee[1]
            if concatee[2] == True :
                numMissing += 1
        outputWriter.write("\n")
        print("strain = " + strainName + ", concatamer length = " + str(totalLength) + ", missing blocks = " + str(numMissing))
    outputWriter.close()
