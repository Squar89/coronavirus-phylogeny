from Bio import Entrez, SeqIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
from matplotlib import pyplot as plt
import os

# create output file where we will write fetched sequences
outputGenomes = open("fetchedGenomes.fa", "w")
outputSpikes = open("fetchedSpikes.fa", "w")

# identifying ourselves to the Entrez database
Entrez.email="jw386401@students.mimuw.edu.pl"
Entrez.tool="testing "

# This function searches and downloads a sequence from Entrez db
def getSequence(database, searchedTerm, name, outputFile, concatWithTerms=[], withWrite=True):
    handle=Entrez.esearch(db=("nuccore" if database == "nucleotide" else "protein"), term=searchedTerm)
    rec=Entrez.read(handle)
    rec_handle=Entrez.efetch(db=("nucleotide" if database == "nucleotide" else "protein"), id=rec["IdList"][0],
                             rettype="fasta")
    sequence=SeqIO.read(rec_handle,"fasta")
    sequence.id=name

    for term in concatWithTerms:
        nextSeq = getSequence(database=database, searchedTerm=term, name="", outputFile=None, withWrite=False)
        sequence.seq += nextSeq.seq

    if withWrite:
        SeqIO.write(sequence, outputFile, "fasta")

    return sequence

# Coronavirus (SARS-CoV-2)
getSequence(database="nucleotide", searchedTerm="MT276331.1", name="SARS-CoV-2", outputFile=outputGenomes)
getSequence(database="protein", searchedTerm="QIS30054.1", name="SARS-CoV-2", outputFile=outputSpikes)

# SARS – another coronavirus back from 2002 (ShanghaiQXC1 strain)
getSequence(database="nucleotide", searchedTerm="AY463059.1", name="SARS", outputFile=outputGenomes)
getSequence(database="protein", searchedTerm="AAR86788.1", name="SARS", outputFile=outputSpikes)

# Bat coronavirus isolate Jiyuan-84
getSequence(database="nucleotide", searchedTerm="KY770860.1", name="Bat-CoV", outputFile=outputGenomes)
getSequence(database="protein", searchedTerm="ARI44809.1", name="Bat-CoV", outputFile=outputSpikes)

# MERS – middle eastern respiratory syndrome coronavirus
getSequence(database="nucleotide", searchedTerm="NC_019843.3", name="MERS", outputFile=outputGenomes)
getSequence(database="protein", searchedTerm="YP_009047204.1", name="MERS", outputFile=outputSpikes)

# Influenza A (H1N1 strain)
getSequence(database="nucleotide", searchedTerm="NC_002023.1", name="Influenza-A", outputFile=outputGenomes,
            concatWithTerms=["NC_002021.1", "NC_002022.1", "NC_002017.1", "NC_002019.1", "NC_002018.1", "NC_002016.1", "NC_002020.1"])

# Hepatitis A
getSequence(database="nucleotide", searchedTerm="NC_001489.1", name="Hepatitis-A", outputFile=outputGenomes)

# Murine hepatitis virus strain A59 (spike)
getSequence(database="protein", searchedTerm="ATN37888.1", name="Murine-Hepatitis", outputFile=outputSpikes)

outputGenomes.close()
outputSpikes.close()

# File path containing genome sequences
fetchedGenomes = "fetchedGenomes.fa"
# File path containing spikes sequences
fetchedSpikes = "fetchedSpikes.fa"

# Generate distance matrix from multiple alignment
def genMatrix(inputFile):
    produceMatrix = ClustalwCommandline("clustalw", infile=inputFile, tree=True, outputtree="dist")
    produceMatrix()

# Generate tree from alignment
def genTree(inputFile, method):
    produceTree = ClustalwCommandline("clustalw", infile=inputFile, tree=True, outputtree="phylip", clustering=method)
    produceTree()

# Perform multiple alignment on genome sequences
alignGenomes = ClustalwCommandline("clustalw", infile=fetchedGenomes, align=True, type="DNA",
                                     pwgapopen=15.00, pwgapext=6.66, gapopen=15.00, gapext=6.66,
                                     maxdiv=30, endgaps=True, transweight=0.50, outorder="aligned",
                                     SEQNOS="ON", SEQNO_RANGE="ON")
alignGenomes()
alignmentGenomes = os.path.splitext(fetchedGenomes)[0] + ".aln" #Get path to the generated alignment file

# Produce distance matrix from multiple alignment results
genMatrix(inputFile=alignmentGenomes)

# Produce tree using NJ algorithm
genTree(inputFile=alignmentGenomes, method="NJ")
# Rename produced file
treeGenomesNJ = os.path.splitext(fetchedGenomes)[0] + ".ph"
os.rename(treeGenomesNJ, os.path.splitext(treeGenomesNJ)[0] + "NJ" + ".ph")
treeGenomesNJ = os.path.splitext(treeGenomesNJ)[0] + "NJ" + ".ph"

#Produce tree using UPGMA algorithm
genTree(inputFile=alignmentGenomes, method="UPGMA")
# Rename produced file
treeGenomesUPGMA = os.path.splitext(fetchedGenomes)[0] + ".ph"
os.rename(treeGenomesUPGMA, os.path.splitext(treeGenomesUPGMA)[0] + "UPGMA" + ".ph")
treeGenomesUPGMA = os.path.splitext(treeGenomesUPGMA)[0] + "UPGMA" + ".ph"

# Perform multiple alignment on spike sequences
alignSpikes = ClustalwCommandline("clustalw", infile=fetchedSpikes, align=True, type="protein", pwgapopen=10.00,
                                  pwgapext=0.10, gapopen=10.00, gapext=0.20, maxdiv=30, endgaps=True, transweight=0.50,
                                  gapdist=4, outorder="aligned", SEQNOS="ON", SEQNO_RANGE="ON")
alignSpikes()
alignmentSpikes = os.path.splitext(fetchedSpikes)[0] + ".aln"

# Produce distance matrix from multiple alignment results
genMatrix(inputFile=alignmentSpikes)

# Produce tree using NJ algorithm
genTree(inputFile=alignmentSpikes, method="NJ")
treeSpikesNJ = os.path.splitext(fetchedSpikes)[0] + ".ph"
os.rename(treeSpikesNJ, os.path.splitext(treeSpikesNJ)[0] + "NJ" + ".ph")
treeSpikesNJ = os.path.splitext(treeSpikesNJ)[0] + "NJ" + ".ph"

# Produce tree using UPGMA algorithm
genTree(inputFile=alignmentSpikes, method="UPGMA")
treeSpikesUPGMA = os.path.splitext(fetchedSpikes)[0] + ".ph"
os.rename(treeSpikesUPGMA, os.path.splitext(treeSpikesUPGMA)[0] + "UPGMA" + ".ph")
treeSpikesUPGMA = os.path.splitext(treeSpikesUPGMA)[0] + "UPGMA" + ".ph"

# Draw produced trees
def drawTreeToFile(tree, outputFilename):
    parsedTree = Phylo.read(tree, "newick")
    Phylo.draw(parsedTree, do_show=False)
    plt.savefig(outputFilename)

drawTreeToFile(treeGenomesNJ, 'genomesTreeNJ.png')
drawTreeToFile(treeGenomesUPGMA, 'genomesTreeUPGMA.png')
drawTreeToFile(treeSpikesNJ, 'spikesTreeNJ.png')
drawTreeToFile(treeSpikesUPGMA, 'spikesTreeUPGMA.png')