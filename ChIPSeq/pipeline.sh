#e li.li@hci.utah.edu
#c kingspeak


NAME=19805X1

# organism Novoalign index
INDEX=/tomato/dev/data/Mouse/Mm10/mm10.standard.nov.illumina.nix

# optical distance
OPTDISTANCE=2500

# adapters for standard compatible Illumina TruSeq
ADAPTF=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
ADATPR=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# application paths
NOVO_APP=/tomato/dev/app/novoalign/4.03.01/novoalign
SAM_APP=/tomato/dev/app/samtools/1.15/samtools
MERGE_APP=/tomato/dev/app/modulesoftware/merge_umi_fastq
DEDUP_APP=/tomato/dev/app/modulesoftware/bam_umi_dedup
BAMWIG_APP=/tomato/dev/app/modulesoftware/bam2wig
MACS_APP=/tomato/dev/app/modulesoftware/macs2
PICARD_APP=/tomato/dev/app/picard/2.23.3/picard.jar

# merge UMI fastq
$MERGE_APP *.fastq.gz

# align
echo "=== aligning"
$NOVO_APP \
-C \
-d $INDEX \
-a $ADAPTF $ADATPR \
-o SAM \
--tune NOVASEQ \
-f *.umi.fastq.gz \
| $SAM_APP fixmate -r -m - $NAME.raw.bam


# sort and index
echo "=== sorting"
$SAM_APP sort -m 4G -@ 12 -o $NAME.sort.bam $NAME.raw.bam \
&& $SAM_APP index -@ 12 $NAME.sort.bam \
&& rm -f $NAME.raw.bam *.umi.fastq.gz


# remove duplicates
echo "=== removing UMI duplicates"
$DEDUP_APP \
--in $NAME.sort.bam \
--out $NAME.bam \
--distance $OPTDISTANCE \
--cpu 12 \
&& rm $NAME.sort.bam*


# check chromosome content
$SAM_APP index -@ 12 $NAME.bam
$SAM_APP idxstats $NAME.bam > $NAME.idxstats.txt


# check insertion sizes
echo "=== insertion sizes"
java -jar $PICARD_APP CollectInsertSizeMetrics \
-I $NAME.bam \
-O $NAME.insertSizeMetrics.txt \
-H $NAME.insertSizes.pdf \
-M 0.5 --VALIDATION_STRINGENCY SILENT --VERBOSITY WARNING 
