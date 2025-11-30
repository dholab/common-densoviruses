#!/usr/bin/env bash
set -euo pipefail

# Generate consensus sequences for BAMs with mapped reads in air-samples/tgs/bam.
# Each of the seven positive samples has exactly one mapped paired-end fragment;
# the consensus is the read-pair sequence oriented to the reference.
# Requirements: samtools (>=1.20) and python3 available on PATH.

BAM_DIR="air-samples/tgs/bam"
CONS_DIR="air-samples/tgs/consensus"

mkdir -p "$CONS_DIR"

echo "Scanning BAMs for mapped reads..."
BAMS=()
for bam in "$BAM_DIR"/*.bam; do
  mapped=$(samtools idxstats "$bam" | awk 'NR==1{print $3}')
  if (( mapped > 0 )); then
    BAMS+=("$bam")
  fi
done

if (( ${#BAMS[@]} == 0 )); then
  echo "No BAMs with mapped reads found. Nothing to do."
  exit 0
fi

for bam in "${BAMS[@]}"; do
  sample=$(basename "$bam" .dedup.sorted.bam)
  # Determine the span covering both mates by parsing CIGAR strings.
  region=$(samtools view -F 4 "$bam" | python3 -c '
import re, sys
min_start = None
max_end = None
chrom = None
pat = re.compile(r"(\d+)([MIDNSHP=X])")
for line in sys.stdin:
    f = line.split("\t")
    chrom = f[2]
    start = int(f[3])
    cigar = f[5]
    span = sum(int(n) for n, op in pat.findall(cigar) if op in "MDN=X")
    end = start + span - 1
    min_start = start if min_start is None else min(min_start, start)
    max_end = end if max_end is None else max(max_end, end)
print(f"{chrom}:{min_start}-{max_end}")
')
  out="$CONS_DIR/${sample}.consensus.fasta"
  echo "Calling consensus for $sample over $region -> $out"
  samtools consensus -r "$region" -f fasta "$bam" | awk -v name="$sample" 'NR==1{print ">"name; next}1' > "$out"
done

echo "Done."
