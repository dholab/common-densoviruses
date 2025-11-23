#!/bin/bash

# Script to index all BAM files for IGV viewing
# This creates .bai index files required for web-based viewing

echo "Indexing BAM files for IGV viewer..."

# Find all BAM files and index them
for bam in air-samples/vsp/bam/*.bam; do
    if [ -f "$bam" ]; then
        bai="${bam}.bai"
        if [ ! -f "$bai" ]; then
            echo "Indexing: $bam"
            samtools index "$bam"
        else
            echo "Index exists: $bai"
        fi
    fi
done

echo "Done! All BAM files are now indexed."