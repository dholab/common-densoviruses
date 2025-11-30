#!/usr/bin/env bash
set -euo pipefail

# Simple end-to-end reproduction of consensus calling, alignment, and tree building.
# Requirements: samtools (>=1.20), MAFFT, IQ-TREE installed and on PATH.
# Reference and BAMs are expected under air-samples/vsp/{ref,bam}.

SAMPLES=(
  elementary-school-004-20250901
  elementary-school-004-20250915
  elementary-school-004-20250922
  elementary-school-004-20250929
  high-school-001-20250901
  high-school-001-20250915
  high-school-001-20250929
)

REF=air-samples/vsp/ref/NC_076998.fasta
BAM_DIR=air-samples/vsp/bam
CONS_DIR=air-samples/vsp/consensus
ALIGN_DIR=air-samples/vsp/align
TREE_DIR=air-samples/vsp/tree

mkdir -p "$CONS_DIR" "$ALIGN_DIR" "$TREE_DIR"

echo "Calling consensus sequences..."
for sample in "${SAMPLES[@]}"; do
  bam="$BAM_DIR/${sample}-pooled.bam"
  out="$CONS_DIR/${sample}-consensus.fasta"
  samtools consensus -f fasta -a -o "$out" "$bam"
done

FA_ALL="$ALIGN_DIR/consensus_plus_ref.fasta"
echo "Combining consensus sequences with reference into $FA_ALL"
: > "$FA_ALL"
{
  echo ">NC_076998_ref"
  awk 'NR>1{print}' "$REF"
  for sample in "${SAMPLES[@]}"; do
    echo ">$sample"
    awk 'NR>1{print}' "$CONS_DIR/${sample}-consensus.fasta"
  done
} >> "$FA_ALL"

ALN="$ALIGN_DIR/consensus_plus_ref.aln.fasta"
echo "Running MAFFT -> $ALN"
# If MAFFT cannot write /dev/stderr (sandboxed), set MAFFT_PROGRESSFILE to a writable path.
: "${MAFFT_PROGRESSFILE:=$ALIGN_DIR/mafft.progress.log}"
MAFFT_PROGRESSFILE="$MAFFT_PROGRESSFILE" mafft --auto --thread 1 "$FA_ALL" > "$ALN"

PREFIX="$TREE_DIR/consensus_plus_ref"
echo "Building IQ-TREE (GTR+G, 1000 UFboot, outgroup NC_076998_ref) -> $PREFIX.*"
iqtree -s "$ALN" -m GTR+G -B 1000 -o NC_076998_ref --prefix "$PREFIX" -T 1

echo "Done."
