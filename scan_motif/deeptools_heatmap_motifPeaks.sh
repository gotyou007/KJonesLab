#!bin/bash


#deeptools_heatmap_motifPeaks

bed=$1
bigwig=$2
bedPrefix=`echo $(basename "$bed" .bed)`
bigwigPrefix=`echo $(basename "$bigwig" .fragment.bw)`
v5=$3

#if v5
if $v5; then
	computeMatrix reference-point \
	  --referencePoint center \
	  -R  $bed \
	  -S $bigwig \
	  -b 4000 -a 4000 \
	  --skipZeros -o matrix1_${bedPrefix}_${bigwigPrefix}.gz \
	  --sortRegions descend \
	  --outFileSortedRegions descending_${bedPrefix}_V5.bed

	plotHeatmap -m matrix1_${bedPrefix}_${bigwigPrefix}.gz -out ${bedPrefix}_${bigwigPrefix}.eps --colorList 'white, blue' --missingDataColor 1 --sortRegions keep --legendLocation none --zMin 0.0 --zMax 2.5

else
	computeMatrix reference-point \
	  --referencePoint center \
	  -R descending_${bedPrefix}_V5.bed \
	  -S $bigwig \
	  -b 4000 -a 4000 \
	  --skipZeros -o matrix1_${bedPrefix}_${bigwigPrefix}.gz \
	  --sortRegions keep
fi
