#e li.li@hci.utah.edu -b
#c kingspeak
#a A6384 # store with the given name on GNomEx

# ref files
ORG=/tomato/dev/data/Human/GRCh38
CHROM=$ORG/chrom.sizes
DB=$ORG/release102
INDEX=$DB/star125
GTF=$DB/Homo_sapiens.GRCh38.102.gtf
REFFLAT=$DB/Homo_sapiens.GRCh38.102.refflat
RIBOINT=$DB/Homo_sapiens.GRCh38.102.rRNA.interval
RSEM_INDEX=$DB/rsem/RSEM

# App Paths
APP=/tomato/dev/app
FASTQC=$APP/FastQC/0.11.5/fastqc
STAR=$APP/STAR/2.7.6a/STAR
FEATCOUNT=$APP/Subread/1.6.3/bin/featureCounts
SAMTOOLS=$APP/samtools/1.10/samtools
BIGWIG=$APP/UCSC/bedGraphToBigWig
PICARD=$APP/picard/2.9.0/picard.jar
RSEM=$APP/rsem/1.3.1/rsem-calculate-expression
CUTADAPT=$APP/modulesoftware/cutadapt
CLUMPIFY=$APP/BBmap/v38.34/clumpify.sh
FASTQSCREEN=$APP/fastq_screen/v0.14.0/fastq_screen

# Concatenate lanes; edit this for your task
R1=`echo *R1_001.fastq.gz`
R2=`echo *R2_001.fastq.gz`

echo $R1
echo $R2

# GET sample prefix
OUT=`echo ${R1%%_*}`

cat $R1 > $OUT.cat_R1.fq.gz
cat $R2 > $OUT.cat_R2.fq.gz

$CLUMPIFY in1=$OUT.cat_R1.fq.gz in2=$OUT.cat_R2.fq.gz out1=$OUT.clump1.fq.gz \
 out2=$OUT.clump2.fq.gz dupedist=12000 dedupe=t optical=t

# CUTADAPT  v 2.8
$CUTADAPT -j $NCPU -O 6 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
 -o $OUT.1.fq -p $OUT.2.fq $OUT.clump1.fq.gz $OUT.clump2.fq.gz

# %rRNA using fastq screen
$FASTQSCREEN --conf /tomato/dev/data/FastQ_Screen_Genomes/eukrRNA/fastq_screen_rRNA.conf \
 --subset 1000000 $OUT.1.fq

rm $OUT.cat_R1.fq.gz $OUT.cat_R2.fq.gz
rm $OUT.clump1.fq.gz $OUT.clump2.fq.gz

# FASTQ
$FASTQC -T $NCPU -f fastq $OUT.1.fq
$FASTQC -T $NCPU -f fastq $OUT.2.fq

# STAR
MEM=`echo "$SMGB" | awk '{print $1 * 1073741824}'`
echo "limitBAMsortRAM=$MEM"
$STAR --genomeDir $INDEX \
--readFilesIn $OUT.1.fq $OUT.2.fq \
--runThreadN $NCPU \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM $MEM \
--outBAMsortingBinsN 100 \
--quantMode TranscriptomeSAM \
--outWigType bedGraph \
--outWigStrand Unstranded

# rename for multiqc ID parsing
mv Aligned.sortedByCoord.out.bam $OUT.bam
mv Log.final.out $OUT.Log.final.out

# Samtools
$SAMTOOLS index $OUT.bam
$SAMTOOLS idxstats $OUT.bam | sort -V > $OUT.idxstats

# RSEM take 'Aligned.toTranscriptome.out' 
$RSEM --paired-end -p $NCPU --alignments --strandedness reverse --no-bam-output \
   Aligned.toTranscriptome.out.bam $RSEM_INDEX $OUT 

# featureCounts -s 2
$FEATCOUNT -T $NCPU -p -s 2 --largestOverlap -a $GTF -o $OUT.counts $OUT.bam
$FEATCOUNT -T $NCPU -p -s 2 --largestOverlap -a $GTF -o $OUT.biotypes -g gene_biotype $OUT.bam

# RnaSeq metrics
/usr/bin/java -Xmx20G -jar $PICARD CollectRnaSeqMetrics  REF_FLAT=$REFFLAT \
STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$RIBOINT \
I=$OUT.bam  O=$OUT.rna_metrics

# bedGraphToBigWig
$BIGWIG Signal.Unique.str1.out.bg $CHROM $OUT.unique.bw
$BIGWIG Signal.UniqueMultiple.str1.out.bg $CHROM $OUT.multiple.bw

LC_COLLATE=C sort -k1,1 -k2,2n Signal.Unique.str1.out.bg > uniq.bg
LC_COLLATE=C sort -k1,1 -k2,2n Signal.UniqueMultiple.str1.out.bg > mult.bg


# keep for Salmon?
rm Aligned.toTranscriptome.out.bam

rm $OUT.1.fq $OUT.2.fq
rm Signal.Unique.str1.out.bg
rm Signal.UniqueMultiple.str1.out.bg
## STAR twopassMode
rm -rf _STARgenome
rm -rf _STARpass1
