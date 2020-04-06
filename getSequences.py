from Bio import Entrez, SeqIO

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