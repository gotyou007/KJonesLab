#!/bin/bash

proj=20828X

for i in {1,4,5};
do
	macs2 callpeak --bdg -q 0.001 -t Alignments/$proj${i}.bam \
		-c Alignments/${proj}6.bam \
	 	-f BAMPE \
		-n $proj$i \
		--outdir macs2
done

for i in {7,8,10,11};
do
	macs2 callpeak --bdg -q 0.001 -t Alignments/$proj${i}.bam \
		-c Alignments/${proj}12.bam \
	 	-f BAMPE \
		-n $proj$i \
		--outdir macs2
done

while IFS="\t" read -r orig new; do echo $orig; done < rename.tsv
	#rename -v "$orig" "$new" *.gz
done < "rename.tsv"