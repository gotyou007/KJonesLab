#!/bin/bash

# For more information about this demo workflow, visit: https://arriba.readthedocs.io/en/latest/workflow/

if [ $# -lt 8 -o $# -gt 9 ]; then
	echo "Usage: $(basename $0) STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv known_fusions.tsv protein_domains.gff3 threads read1.fastq.gz [read2.fastq.gz]" 1>&2
	exit 1
fi

# tell bash to be verbose and to abort on error
set -o pipefail
set -x -e -u

# get arguments
STAR_INDEX_DIR="$1"
ANNOTATION_GTF="$2"
ASSEMBLY_FA="$3"
BLACKLIST_TSV="$4"
KNOWN_FUSIONS_TSV="$5"
TAGS_TSV="$KNOWN_FUSIONS_TSV" # different files can be used for filtering and tagging, but the provided one can be used for both
PROTEIN_DOMAINS_GFF3="$6"
THREADS="$7"
READ1="$8"
READ2="${9-}"

# find installation directory of arriba
BASE_DIR=$(dirname "$0")

# align FastQ files (STAR >=2.7.10a recommended)
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn "$READ1" "$READ2" --readFilesCommand gzcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |

tee Aligned.out.bam |

# call arriba
"$BASE_DIR/arriba" \
	-x /dev/stdin \
	-o fusions.tsv -O fusions.discarded.tsv \
	-a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" -k "$KNOWN_FUSIONS_TSV" -t "$TAGS_TSV" -p "$PROTEIN_DOMAINS_GFF3" \
#	-d structural_variants_from_WGS.tsv

# sorting and indexing is only required for visualization
if [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	samtools sort -@ "$THREADS" -m $((40000/THREADS))M -T tmp -O bam Aligned.out.bam > Aligned.sortedByCoord.out.bam
	rm -f Aligned.out.bam
	samtools index Aligned.sortedByCoord.out.bam
else
	echo "samtools >= 1.0 required for sorting of alignments" 1>&2
fi

