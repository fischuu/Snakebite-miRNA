# tRNA

Get the tRNA mature database from here:
https://gtrnadb.ucsc.edu/index.html

For example, for bos taurus, take the mature fasta file from here:
https://gtrnadb.ucsc.edu/genomes/eukaryota/Btaur9/bosTau9-mature-tRNAs.fa

# PhiX
wget https://webdata.illumina.com/downloads/productfiles/igenomes/phix/PhiX_Illumina_RTA.tar.gz
tar -xzvf PhiX_Illumina_RTA.tar.gz
cd PhiX/Illumina/RTA/Sequence/WholeGenomeFasta
cp genome.fa $PROJECTFOLDER/references/phix.fa

# Mirbase
wget https://www.mirbase.org/download/mature.fa
wget https://www.mirbase.org/download/hairpin.fa