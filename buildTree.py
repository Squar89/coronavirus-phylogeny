from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
import os

# File containing genome sequences
fetchedGenomes = "fetchedGenomes.fa"

# Perform multiple alignment on genome sequences
alignGenomes = ClustalwCommandline("clustalw", infile=fetchedGenomes, align=True, type="DNA",
                                     pwgapopen=15.00, pwgapext=6.66, gapopen=15.00, gapext=6.66,
                                     maxdiv=30, endgaps=True, transweight=0.50, outorder="aligned",
                                     SEQNOS="ON", SEQNO_RANGE="ON")
print(alignGenomes)
alignGenomes()
alignmentGenomes = os.path.splitext(fetchedGenomes)[0] + ".aln"

# Produce distance matrix from multiple alignment results
produceMatrixGenomes = ClustalwCommandline("clustalw", infile=alignmentGenomes, tree=True, outputtree="dist")
print(produceMatrixGenomes)
produceMatrixGenomes()

# Produce tree using NJ algorithm
treeGenomesNJ = os.path.splitext(fetchedGenomes)[0] + ".ph"
produceTreeGenomesNJ = ClustalwCommandline("clustalw", infile=alignmentGenomes, outfile=treeGenomesNJ, tree=True,
                                           outputtree="phylip", clustering="NJ")
print(produceTreeGenomesNJ)
produceTreeGenomesNJ()
os.rename(treeGenomesNJ, os.path.splitext(treeGenomesNJ)[0] + "NJ" + ".ph")
treeGenomesNJ = os.path.splitext(treeGenomesNJ)[0] + "NJ" + ".ph"

#Produce tree using UPGMA algorithm
treeGenomesUPGMA = os.path.splitext(fetchedGenomes)[0] + ".ph"
produceTreeGenomesUPGMA = ClustalwCommandline("clustalw", infile=alignmentGenomes, outfile=treeGenomesUPGMA, tree=True,
                                           outputtree="phylip", clustering="UPGMA")
print(produceTreeGenomesUPGMA)
produceTreeGenomesUPGMA()
os.rename(treeGenomesUPGMA, os.path.splitext(treeGenomesUPGMA)[0] + "UPGMA" + ".ph")
treeGenomesUPGMA = os.path.splitext(treeGenomesUPGMA)[0] + "UPGMA" + ".ph"

# Repeat the process for spike proteins
fetchedSpikes = "fetchedSpikes.fa"

# Perform multiple alignment on spike sequences
alignSpikes = ClustalwCommandline("clustalw", infile=fetchedSpikes, align=True, type="protein", pwgapopen=10.00,
                                  pwgapext=0.10, gapopen=10.00, gapext=0.20, maxdiv=30, endgaps=True, transweight=0.50,
                                  gapdist=4, outorder="aligned", SEQNOS="ON", SEQNO_RANGE="ON")

print(alignSpikes)
alignSpikes()
alignmentSpikes = os.path.splitext(fetchedSpikes)[0] + ".aln"

# Produce distance matrix from multiple alignment results
produceMatrixSpikes = ClustalwCommandline("clustalw", infile=alignmentSpikes, tree=True, outputtree="dist")
print(produceMatrixSpikes)
produceMatrixSpikes()

# Produce tree using NJ algorithm
treeSpikesNJ = os.path.splitext(fetchedSpikes)[0] + ".ph"
produceTreeSpikesNJ = ClustalwCommandline("clustalw", infile=alignmentSpikes, tree=True, outputtree="phylip",
                                          clustering="NJ")
print(produceTreeSpikesNJ)
produceTreeSpikesNJ()
os.rename(treeSpikesNJ, os.path.splitext(treeSpikesNJ)[0] + "NJ" + ".ph")
treeSpikesNJ = os.path.splitext(treeSpikesNJ)[0] + "NJ" + ".ph"

# Produce tree using UPGMA algorithm
treeSpikesUPGMA = os.path.splitext(fetchedSpikes)[0] + ".ph"
produceTreeSpikesUPGMA = ClustalwCommandline("clustalw", infile=alignmentSpikes, tree=True, outputtree="phylip",
                                             clustering="UPGMA")
print(produceTreeSpikesUPGMA)
produceTreeSpikesUPGMA()
os.rename(treeSpikesUPGMA, os.path.splitext(treeSpikesUPGMA)[0] + "UPGMA" + ".ph")
treeSpikesUPGMA = os.path.splitext(treeSpikesUPGMA)[0] + "UPGMA" + ".ph"

# Draw produced trees
phyloTreeGenomesNJ = Phylo.read(treeGenomesNJ, "newick")
Phylo.draw_ascii(phyloTreeGenomesNJ)

phyloTreeGenomesUPGMA = Phylo.read(treeGenomesUPGMA, "newick")
Phylo.draw_ascii(phyloTreeGenomesUPGMA)

phyloTreeSpikesNJ = Phylo.read(treeSpikesNJ, "newick")
Phylo.draw_ascii(phyloTreeSpikesNJ)

phyloTreeSpikesUPGMA = Phylo.read(treeSpikesUPGMA, "newick")
Phylo.draw_ascii(phyloTreeSpikesUPGMA)