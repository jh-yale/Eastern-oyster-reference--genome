#!/bin/bash
# RNA-seq guided gene prediction on the removed sequences, reusing the
# Augustus species model from the original genome annotation so the
# comparison is apples-to-apples

set -euo pipefail

module load GCC/13.3.0 STAR SAMtools

READS1=$1
READS2=$2
SPECIES=$3
THREADS=16

WORKDIR=manual_removal
QUERY=$WORKDIR/removed_sequences_nonN.fa
ANNOTDIR=$WORKDIR/annotation
mkdir -p $ANNOTDIR

mkdir -p $ANNOTDIR/star_index
STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir $ANNOTDIR/star_index \
     --genomeFastaFiles $QUERY \
     --genomeSAindexNbases 8

STAR --runThreadN $THREADS \
     --genomeDir $ANNOTDIR/star_index \
     --readFilesIn $READS1 $READS2 \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 4000000000 \
     --outFileNamePrefix $ANNOTDIR/rnaseq_

samtools index $ANNOTDIR/rnaseq_Aligned.sortedByCoord.out.bam

[ -f $ANNOTDIR/braker3.sif ] || apptainer pull $ANNOTDIR/braker3.sif docker://teambraker/braker3:latest

apptainer exec \
  --bind $WORKDIR:$WORKDIR \
  --bind $AUGUSTUS_CONFIG_PATH:/config \
  $ANNOTDIR/braker3.sif \
  bash -c "
    export AUGUSTUS_CONFIG_PATH=/config
    export GENEMARK_PATH=/opt/ETP/bin
    export AUGUSTUS_BIN_PATH=/opt/Augustus/bin
    export AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts
    export PYTHON3_PATH=/opt/conda/bin

    braker.pl --genome=$QUERY --bam=$ANNOTDIR/rnaseq_Aligned.sortedByCoord.out.bam \
      --softmasking --threads=$THREADS --species=$SPECIES \
      --workingdir=$ANNOTDIR/braker_output --gff3 --useexisting
  "
