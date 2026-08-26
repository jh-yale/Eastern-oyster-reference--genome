#!/bin/bash
# check whether predicted genes on the removed sequences are represented
# in the final assembly -- first at the protein level, then (for anything
# that doesn't match a protein) at the raw genomic DNA level

set -euo pipefail

module load GCC/13.3.0 BLAST+ SAMtools

FINAL_PROTEOME=$1
FINAL_ASSEMBLY=$2
THREADS=8

WORKDIR=manual_removal
ANNOTDIR=$WORKDIR/annotation
GFF=$ANNOTDIR/braker_output/braker.gff3
PROTEINS=$ANNOTDIR/braker_output/braker.aa
COMPAREDIR=$WORKDIR/comparison
mkdir -p $COMPAREDIR

# --- protein-level check ---
makeblastdb -in "$FINAL_PROTEOME" -dbtype prot -out $COMPAREDIR/final_proteins_db

blastp -query $PROTEINS -db $COMPAREDIR/final_proteins_db \
    -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
    -evalue 1e-5 -num_threads $THREADS -max_target_seqs 1 \
    -out $COMPAREDIR/proteins_vs_final_blastp.txt

grep "^>" $PROTEINS | sed 's/>//' | sort > $COMPAREDIR/all_ids.txt
cut -f1 $COMPAREDIR/proteins_vs_final_blastp.txt | sort -u > $COMPAREDIR/matched_ids.txt
comm -23 $COMPAREDIR/all_ids.txt $COMPAREDIR/matched_ids.txt > $COMPAREDIR/nohit_ids.txt

echo "total predictions: $(wc -l < $COMPAREDIR/all_ids.txt)"
echo "matched a final-assembly protein: $(wc -l < $COMPAREDIR/matched_ids.txt)"
echo "no protein match: $(wc -l < $COMPAREDIR/nohit_ids.txt)"

[ -s $COMPAREDIR/nohit_ids.txt ] || exit 0

# --- for the no-hit set, pull the unspliced gene region and blast that
# against the whole final assembly instead of just its annotated proteins ---
> $COMPAREDIR/nohit_gene_coords.tsv
while read tid; do
    gid=$(echo "$tid" | sed -E 's/\.t[0-9]+$//')
    line=$(awk -F'\t' -v g="ID=${gid};" '$3=="gene" && index($9,g)' "$GFF" | head -1)
    [ -z "$line" ] && continue
    echo -e "$tid\t$gid\t$(echo "$line" | cut -f1)\t$(echo "$line" | cut -f4)\t$(echo "$line" | cut -f5)"
done < $COMPAREDIR/nohit_ids.txt > $COMPAREDIR/nohit_gene_coords.tsv

[ -f $WORKDIR/removed_sequences_nonN.fa.fai ] || samtools faidx $WORKDIR/removed_sequences_nonN.fa

> $COMPAREDIR/nohit_genomic_regions.fa
while IFS=$'\t' read -r tid gid src start end; do
    samtools faidx $WORKDIR/removed_sequences_nonN.fa "${src}:${start}-${end}" 2>/dev/null |
        sed "1s/.*/>${tid}/" \
        >> $COMPAREDIR/nohit_genomic_regions.fa
done < $COMPAREDIR/nohit_gene_coords.tsv

makeblastdb -in "$FINAL_ASSEMBLY" -dbtype nucl -out $COMPAREDIR/final_genome_db

blastn -query $COMPAREDIR/nohit_genomic_regions.fa -db $COMPAREDIR/final_genome_db \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore" \
    -evalue 1e-10 -num_threads $THREADS -max_target_seqs 3 \
    -out $COMPAREDIR/nohit_genomic_vs_final_blastn.txt

sort -k1,1 -k8,8rn $COMPAREDIR/nohit_genomic_vs_final_blastn.txt |
    awk '!seen[$1]++' > $COMPAREDIR/nohit_genomic_besthit.txt

TOTAL=$(grep -c "^>" $COMPAREDIR/nohit_genomic_regions.fa)
WITH_DNA=$(cut -f1 $COMPAREDIR/nohit_genomic_besthit.txt | sort -u | wc -l)

echo ""
echo "protein-unmatched predictions checked at the DNA level: $TOTAL"
echo "genomic sequence found in final assembly: $WITH_DNA"
echo "no match anywhere: $((TOTAL - WITH_DNA))"
awk '{i+=$3; c+=$6; n++} END{if(n>0) printf "mean identity %.1f%%, mean coverage %.1f%%\n", i/n, c/n}' \
    $COMPAREDIR/nohit_genomic_besthit.txt

comm -23 \
    <(grep "^>" $COMPAREDIR/nohit_genomic_regions.fa | sed 's/>//' | sort) \
    <(cut -f1 $COMPAREDIR/nohit_genomic_besthit.txt | sort -u)
