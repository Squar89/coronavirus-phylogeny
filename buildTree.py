from Bio.Align.Applications import ClustalwCommandline

in_file = "fetchedGenomes.fa"
clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
print(clustalw_cline)
clustalw_cline()