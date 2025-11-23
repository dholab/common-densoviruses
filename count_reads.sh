#!/bin/bash
for bam in docs/data/bam/*.bam; do
    count=$(samtools view -c "$bam" 2>/dev/null)
    if [ "$count" -gt 0 ]; then
        name=$(basename "$bam" .bam)
        echo "$name: $count"
    fi
done | sort -t: -k2 -rn