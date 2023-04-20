#!/bin/bash

docker run --rm -v `pwd`:/references uhrigs/arriba:2.4.0 download_references.sh hg38+RefSeq


./run_arriba.sh references/star125 references/RefSeq_hg38.gtf references/hg38.fa database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz database/protein_domains_hg38_GRCh38_v2.4.0.gff3 12 /Users/li/bio_info/A7968/A7968/pdxRun/read1.fastq.gz /Users/li/bio_info/A7968/A7968/pdxRun/read2.fastq.gz


./draw_fusions.R \
    --fusions=fusions.tsv \
    --output=fusions.pdf \
    --annotation=references/RefSeq_hg38.gtf \
    --cytobands=database/cytobands_hg38_GRCh38_v2.4.0.tsv \
    --proteinDomains=database/protein_domains_hg38_GRCh38_v2.4.0.gff3

/Users/li/detect-fusion/arriba_v2.4.0/run_arriba.sh /Users/li/detect-fusion/arriba_v2.4.0/references/star125 \
  /Users/li/detect-fusion/arriba_v2.4.0/references/RefSeq_hg38.gtf \
  /Users/li/detect-fusion/arriba_v2.4.0/references/hg38.fa \
  /Users/li/detect-fusion/arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz \
  /Users/li/detect-fusion/arriba_v2.4.0/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz \
  /Users/li/detect-fusion/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 12 \
  20643X9_230324_A00421_0538_AH2TC5DSX7_S9_L001_R1_001.fastq.gz \
  20643X9_230324_A00421_0538_AH2TC5DSX7_S9_L001_R2_001.fastq.gz

/Users/li/detect-fusion/arriba_v2.4.0/draw_fusions.R \
    --fusions=fusions.tsv \
    --output=fusions.pdf \
    --annotation=/Users/li/detect-fusion/arriba_v2.4.0/references/RefSeq_hg38.gtf \
    --cytobands=/Users/li/detect-fusion/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv \
    --proteinDomains=/Users/li/detect-fusion/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3