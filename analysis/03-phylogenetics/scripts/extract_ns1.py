#!/usr/bin/env python3
"""
Extract NS1 amino acid sequences from GenBank file.
Searches for product names containing NS1, non-structural protein 1,
or nonstructural protein 1 (case-insensitive).
"""

from Bio import SeqIO
import re
import sys

def is_ns1_product(product_name):
    """Check if product name indicates NS1 protein."""
    if not product_name:
        return False
    product_lower = product_name.lower()
    # Match various NS1 naming conventions:
    # - NS1, NS-1, NS 1
    # - non-structural protein 1, nonstructural protein 1
    # - nonstructural protein NS1, non-structural protein NS-1, etc.
    # - helicase NS1, NS1 protein, putative NS1, etc.
    patterns = [
        r'\bns[- ]?1\b',                      # NS1, NS-1, NS 1
        r'non-structural protein 1\b',        # non-structural protein 1
        r'nonstructural protein 1\b',         # nonstructural protein 1
        r'non-structural protein ns[- ]?1',   # non-structural protein NS1/NS-1
        r'nonstructural protein ns[- ]?1',    # nonstructural protein NS1/NS-1
        r'helicase ns[- ]?1',                 # helicase NS1
    ]
    for pattern in patterns:
        if re.search(pattern, product_lower):
            return True
    return False

def sanitize_header(text):
    """Replace whitespaces with underscores."""
    if text:
        return re.sub(r'\s+', '_', text.strip())
    return ""

def extract_ns1_sequences(genbank_file, output_file):
    """Extract NS1 amino acid sequences from GenBank file."""

    ns1_records = []

    for record in SeqIO.parse(genbank_file, "genbank"):
        nucleotide_accession = record.id

        # Get host from source feature
        host = ""
        for feature in record.features:
            if feature.type == "source":
                if "host" in feature.qualifiers:
                    host = feature.qualifiers["host"][0]
                break

        # Search for CDS features with NS1 product
        for feature in record.features:
            if feature.type == "CDS":
                # Get product name
                product = feature.qualifiers.get("product", [""])[0]
                gene = feature.qualifiers.get("gene", [""])[0]

                # Check if this is NS1
                if is_ns1_product(product) or is_ns1_product(gene):
                    # Get protein_id
                    protein_id = feature.qualifiers.get("protein_id", [""])[0]

                    # Get translation (amino acid sequence)
                    translation = feature.qualifiers.get("translation", [""])[0]

                    if translation:
                        # Create header: [nucleotide accession]_[protein accession]_[host]
                        header_parts = [
                            sanitize_header(nucleotide_accession),
                            sanitize_header(protein_id),
                            sanitize_header(host)
                        ]
                        header = "_".join(part for part in header_parts if part)

                        ns1_records.append((header, translation, product))

    # Write to FASTA file
    with open(output_file, 'w') as f:
        for header, seq, product in ns1_records:
            f.write(f">{header}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    return len(ns1_records)

if __name__ == "__main__":
    genbank_file = "phylogenetics/known densoviruses NCBI genbank 2025-11-19.gb"
    output_file = "phylogenetics/ns1_sequences.fasta"

    count = extract_ns1_sequences(genbank_file, output_file)
    print(f"Extracted {count} NS1 sequences to {output_file}")
