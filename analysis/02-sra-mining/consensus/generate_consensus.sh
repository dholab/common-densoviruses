#!/usr/bin/env bash
set -euo pipefail

# Generate consensus sequences from Logan BAM files and identify those with >=80% coverage

REF_LEN=4678
MIN_COV=0.80
MIN_BASES=$(echo "$REF_LEN * $MIN_COV" | bc | cut -d. -f1)

BAM_DIR="docs/02-sra-mining/data/bam"
CONS_DIR="analysis/02-sra-mining/consensus"

echo "Reference length: $REF_LEN bp"
echo "Minimum coverage threshold: 80% = $MIN_BASES bases"
echo ""

# Clear previous results
> "$CONS_DIR/coverage_stats.txt"
> "$CONS_DIR/genomes.txt"

echo "Generating consensus sequences and calculating coverage..."
for bam in $BAM_DIR/*.dedup.sorted.bam; do
  sample=$(basename "$bam" .dedup.sorted.bam)
  out="$CONS_DIR/${sample}.consensus.fasta"

  # Generate consensus with samtools consensus
  samtools consensus -f fasta -a -o "$out" "$bam" 2>/dev/null

  # Count non-N bases in consensus (excluding header line)
  covered=$(awk 'NR>1' "$out" | tr -d '\n' | sed 's/[nN*]//g' | wc -c | tr -d ' ')
  pct=$(echo "scale=2; $covered * 100 / $REF_LEN" | bc)

  echo "$sample: $covered bases covered (${pct}%)" >> "$CONS_DIR/coverage_stats.txt"

  # Check if coverage meets threshold
  if [ "$covered" -ge "$MIN_BASES" ]; then
    echo "$sample" >> "$CONS_DIR/genomes.txt"
  fi
done

echo ""
echo "Coverage stats saved to $CONS_DIR/coverage_stats.txt"
echo ""
echo "Samples with >=80% coverage:"
cat "$CONS_DIR/genomes.txt"
echo ""
echo "Total samples with >=80% coverage: $(wc -l < "$CONS_DIR/genomes.txt")"
