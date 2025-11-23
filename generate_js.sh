#!/bin/bash
echo "        const samples = ["
for bam in docs/data/bam/*.bam; do
    count=$(samtools view -c "$bam" 2>/dev/null)
    name=$(basename "$bam" .bam)
    echo "$count:$name"
done | sort -t: -k1 -rn | while IFS=: read count name; do
    echo "            { name: \"$name\", reads: $count },"
done | sed '$ s/,$//'
echo "        ];"