# Script to pick fasta entries of interest out of the full DB dump of GISAID DB
# The obtained fasta files are merged and outputted as a new fasta file 
# The output fasta file is subsequently used for the indexing of kraken2 DB

#!/bin/bash


targetDir=~/scratch/gisaid
k2FASTAfile=$targetDir/relevantGISAID_k2.fa
kallistoFASTAfile=$targetDir/relevantGISAID_kallisto.fa
taxonomyFile=./majorCovidDB/taxonomy/names.dmp


echo ">wt-wuhan|kraken:taxid|1000" > $k2FASTAfile
cat ../covidRefSequences/wuhan.fa | grep -v '>' >> $k2FASTAfile
echo >> $k2FASTAfile

echo ">wt-wuhan|kraken:taxid|1000" > $kallistoFASTAfile
cat ../covidRefSequences/wuhan.fa | grep -v '>' >> $kallistoFASTAfile
echo >> $kallistoFASTAfile


for gisaidID in `cat $taxonomyFile | grep "EPI_ISL" | awk -F '|' '{print $2 }'`; do
	metadataRow=`cat $targetDir/metadata.tsv | grep -w $gisaidID`	
	fastaHeader=`echo "$metadataRow" | awk -F $'\t' '{ print $1 }'`
	pangoLineage=`echo "$metadataRow" | awk -F $'\t' '{ print $12 }'`
	echo $gisaidID "$fastaHeader" $pangoLineage
	
	grep -A 1000 -m 1 "$fastaHeader" $targetDir/sequences.fasta > $targetDir/excerpt.fa
	head -n 1 $targetDir/excerpt.fa
	csplit -kzq --prefix $targetDir/xx $targetDir/excerpt.fa "/>/" "{*}"
	
	krakenTaxid=`cat $taxonomyFile | grep $gisaidID | awk '{ print $1 }'`
	echo ">$gisaidID|kraken:taxid|$krakenTaxid" >> $k2FASTAfile
	cat $targetDir/xx00 | grep -v '>' >> $k2FASTAfile
	echo >> $k2FASTAfile
	
	echo ">$pangoLineage" >> $kallistoFASTAfile
	cat $targetDir/xx00 | grep -v '>' >> $kallistoFASTAfile
	echo >> $kallistoFASTAfile
	rm $targetDir/xx0*
done

# Index the allCovidDB for kraken2
kraken2-build --db allCovidDB --add-to-library $k2FASTAfile 
kraken2-build --build --db allCovidDB
bracken-build -d allCovidDB -t 10

# Index the majorCovidDB for kraken2
kraken2-build --db majorCovidDB --add-to-library $k2FASTAfile 
kraken2-build --build --db majorCovidDB
bracken-build -d majorCovidDB -t 10

# Generate kallisto index
kallisto index --index variants.kalIdx $kallistoFASTAfile --make-unique


