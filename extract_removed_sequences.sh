#!/bin/bash
# identify sequence present in the pre-curation assembly but missing from
# the final assembly (i.e. removed during manual haplotig cleanup)

set -euo pipefail

module load GCC/13.3.0 minimap2 SAMtools BEDTools

PRE=$1
FINAL=$2
THREADS=16

WORKDIR=manual_removal
mkdir -p $WORKDIR

ln -sf "$(readlink -f "$PRE")" $WORKDIR/pre_assembly.fa
samtools faidx $WORKDIR/pre_assembly.fa
cut -f1,2 $WORKDIR/pre_assembly.fa.fai | sort -k1,1 > $WORKDIR/pre_assembly.sizes

# align pre-curation assembly onto the final assembly
minimap2 -ax asm5 -t $THREADS "$FINAL" $WORKDIR/pre_assembly.fa 2> $WORKDIR/minimap2.log |
    samtools sort -@ $THREADS -o $WORKDIR/pre_on_final.bam
samtools index $WORKDIR/pre_on_final.bam

bedtools bamtobed -i $WORKDIR/pre_on_final.bam |
    sort -k1,1 -k2,2n |
    bedtools merge -i stdin \
    > $WORKDIR/pre_aligned.bed

# complement of aligned regions = removed sequence
bedtools complement -i $WORKDIR/pre_aligned.bed -g $WORKDIR/pre_assembly.sizes \
    > $WORKDIR/removed_regions_raw.bed

awk '($3-$2) >= 200' $WORKDIR/removed_regions_raw.bed > $WORKDIR/removed_regions.bed

bedtools getfasta -fi $WORKDIR/pre_assembly.fa -bed $WORKDIR/removed_regions.bed \
    -fo $WORKDIR/removed_sequences.fa

# drop pure-N fragments (assembly gap placeholders, not real sequence)
awk '
/^>/ { if (h != "" && s !~ /^N+$/) { print h; print s }
       h=$0; s=""; next }
{ s = s $0 }
END { if (h != "" && s !~ /^N+$/) { print h; print s } }
' $WORKDIR/removed_sequences.fa > $WORKDIR/removed_sequences_nonN.fa

# sanity check -- each removed fragment should map back to its own
# labeled position in the pre-curation assembly at ~100% identity/coverage
minimap2 -x asm5 -t $THREADS $WORKDIR/pre_assembly.fa $WORKDIR/removed_sequences_nonN.fa 2>/dev/null \
    > $WORKDIR/removed_selfcheck.paf
awk '{print $1, $2, ($4-$3)/$2*100"% coverage"}' $WORKDIR/removed_selfcheck.paf

echo "fragments (all): $(grep -c '^>' $WORKDIR/removed_sequences.fa)"
echo "fragments (non-N): $(grep -c '^>' $WORKDIR/removed_sequences_nonN.fa)"
echo "total bp: $(awk '/^>/{next}{n+=length($0)}END{print n}' $WORKDIR/removed_sequences_nonN.fa)"
