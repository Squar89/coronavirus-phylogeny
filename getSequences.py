from Bio import Entrez
from Bio import SeqIO

# create output file where we will write fetched sequences
outputGenomes = open("fetchedGenomes.fa", "w")
outputSpikes = open("fetchedSpikes.fa", "w")

# identifying ourselves to the Entrez database
Entrez.email="jw386401@students.mimuw.edu.pl"
Entrez.tool="testing "

# searching for coronavirus sequence
handle=Entrez.esearch(db="nuccore",term="2019-nCov")
rec=Entrez.read(handle)

# now we have the search results in a dictionary, we can take the list of sequence IDs
rec["IdList"]

#we can fetch the first one - the whole genome sequence
rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][2],rettype="fasta")

# now we have the handle, we need to read the fasta file from it 
ncov19=SeqIO.read(rec_handle,"fasta")

# write found sequence to output file
ncov19.id = "SARS-CoV-2"
SeqIO.write(ncov19, outputGenomes, "fasta")

# now we will do the same for the spike protein
handle=Entrez.esearch(db="protein",term="2019-nCov spike protein")
spike_rec=Entrez.read(handle)
spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][0],rettype="fasta")
spike_ncov19=SeqIO.read(spike_handle,"fasta")
spike_ncov19.id = "SARS-CoV-2"
SeqIO.write(spike_ncov19, outputSpikes, "fasta")

# repeat the same process for other sequences
# SARS – another coronavirus back from 2002 (ShanghaiQXC1 strain)
handle=Entrez.esearch(db="nuccore",term="AY463059.1")
rec=Entrez.read(handle)
rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][0],rettype="fasta")
sars1=SeqIO.read(rec_handle,"fasta")
sars1.id = "SARS"
SeqIO.write(sars1, outputGenomes, "fasta")

handle=Entrez.esearch(db="protein",term="AAR86788.1")
spike_rec=Entrez.read(handle)
spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][0],rettype="fasta")
spike_sars1=SeqIO.read(spike_handle,"fasta")
spike_sars1.id = "SARS"
SeqIO.write(spike_sars1, outputSpikes, "fasta")

# Bat coronavirus isolate Jiyuan-84, complete genome
handle=Entrez.esearch(db="nuccore",term="KY770860.1")
rec=Entrez.read(handle)
rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][0],rettype="fasta")
batCov=SeqIO.read(rec_handle,"fasta")
batCov.id = "Bat-CoV"
SeqIO.write(batCov, outputGenomes, "fasta")

handle=Entrez.esearch(db="protein",term="ARI44809.1")
spike_rec=Entrez.read(handle)
spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][0],rettype="fasta")
spike_batCov=SeqIO.read(spike_handle,"fasta")
spike_batCov.id = "Bat-CoV"
SeqIO.write(spike_batCov, outputSpikes, "fasta")

# MERS – middle eastern respiratory syndrome coronavirus
handle=Entrez.esearch(db="nuccore",term="NC_019843.3")
rec=Entrez.read(handle)
rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][0],rettype="fasta")
mers=SeqIO.read(rec_handle,"fasta")
mers.id = "MERS"
SeqIO.write(mers, outputGenomes, "fasta")

handle=Entrez.esearch(db="protein",term="YP_009047204.1")
spike_rec=Entrez.read(handle)
spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][0],rettype="fasta")
spike_mers=SeqIO.read(spike_handle,"fasta")
spike_mers.id = "MERS"
SeqIO.write(spike_mers, outputSpikes, "fasta")

# Influenza A (H1N1 strain)
handle=Entrez.esearch(db="nuccore",term="NC_002023.1")#First segment
rec=Entrez.read(handle)
rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][0],rettype="fasta")
influenza=SeqIO.read(rec_handle,"fasta")

segmentsIdList = ["NC_002021.1", "NC_002022.1", "NC_002017.1", "NC_002019.1", "NC_002018.1", "NC_002016.1", "NC_002020.1"]
for id in segmentsIdList:
    handle=Entrez.esearch(db="nuccore",term=id)
    rec=Entrez.read(handle)
    rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][0],rettype="fasta")
    influenzaNextSegment=SeqIO.read(rec_handle,"fasta")
    influenza.seq += influenzaNextSegment.seq

influenza.id = "Influenza-A"
influenza.description = influenza.description.replace(" segment 1, complete sequence", ", complete genome")
# Write all segments to the file
SeqIO.write(influenza, outputGenomes, "fasta")

# Hepatitis A virus
handle=Entrez.esearch(db="nuccore",term="NC_001489.1")
rec=Entrez.read(handle)
rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][0],rettype="fasta")
hepatitis=SeqIO.read(rec_handle,"fasta")
hepatitis.id = "Hepatitis-A"
SeqIO.write(hepatitis, outputGenomes, "fasta")

# Murine hepatitis virus strain A59 (spike)
handle=Entrez.esearch(db="protein",term="ATN37888.1")
spike_rec=Entrez.read(handle)
spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][0],rettype="fasta")
spike_hepatitis=SeqIO.read(spike_handle,"fasta")
spike_hepatitis.id = "Murine-Hepatitis"
SeqIO.write(spike_hepatitis, outputSpikes, "fasta")