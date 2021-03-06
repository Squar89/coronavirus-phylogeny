from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from matplotlib import pyplot as plt
import os

# File containing genome sequences
fetchedGenomes = "fetchedGenomes.fa"
# File containing spikes sequences
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