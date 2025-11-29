#!/usr/bin/env python3
"""
Extract NS1 amino acid sequences from VSP consensus files.
Uses GFF annotation to determine NS1 CDS coordinates and translates to amino acids.
"""

import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_ns1_aa_sequences():
    """Extract NS1 amino acid sequences from all consensus files."""

    # NS1 CDS coordinates from GFF file (1-based, converting to 0-based for Python)
    # CDS 238-301 and CDS 373-2015
    ns1_cds1_start = 238 - 1  # 237 (0-based)
    ns1_cds1_end = 301       # 301 (0-based, exclusive)
    ns1_cds2_start = 373 - 1  # 372 (0-based)
    ns1_cds2_end = 2015      # 2015 (0-based, exclusive)

    # Find all consensus files
    consensus_files = glob.glob("air-samples/vsp/consensus/*.fasta")
    consensus_files.sort()

    # List to store all NS1 amino acid sequences
    ns1_aa_records = []

    for fasta_file in consensus_files:
        print(f"Processing {fasta_file}")

        # Extract sample name from filename
        filename = os.path.basename(fasta_file)
        sample_name = filename.replace("-consensus.fasta", "")

        # Read the consensus sequence
        try:
            record = next(SeqIO.parse(fasta_file, "fasta"))
            dna_seq = str(record.seq)

            # Extract NS1 CDS regions
            ns1_cds1 = dna_seq[ns1_cds1_start:ns1_cds1_end]
            ns1_cds2 = dna_seq[ns1_cds2_start:ns1_cds2_end]

            # Concatenate the two CDS regions
            ns1_complete_cds = ns1_cds1 + ns1_cds2

            print(f"  NS1 CDS1 length: {len(ns1_cds1)} bp")
            print(f"  NS1 CDS2 length: {len(ns1_cds2)} bp")
            print(f"  Complete NS1 CDS length: {len(ns1_complete_cds)} bp")

            # Check if sequence length is divisible by 3
            if len(ns1_complete_cds) % 3 != 0:
                print(f"  Warning: NS1 CDS length ({len(ns1_complete_cds)}) not divisible by 3 for {sample_name}")
                # Trim to the largest multiple of 3
                ns1_complete_cds = ns1_complete_cds[:len(ns1_complete_cds) - (len(ns1_complete_cds) % 3)]
                print(f"  Trimmed to {len(ns1_complete_cds)} bp")

            # Convert to Biopython Seq object and translate
            ns1_dna_seq = Seq(ns1_complete_cds)

            # Translate to amino acids, handling ambiguous nucleotides properly
            try:
                # First, let's handle ambiguous nucleotides by translating codon by codon
                aa_sequence = ""
                for i in range(0, len(ns1_complete_cds), 3):
                    codon = ns1_complete_cds[i:i+3]
                    if len(codon) == 3:
                        # If codon contains any ambiguous nucleotides, use 'X'
                        if 'N' in codon or len(codon) < 3:
                            aa_sequence += 'X'
                        else:
                            try:
                                # Translate the clean codon
                                aa = Seq(codon).translate()
                                aa_sequence += str(aa)
                            except:
                                aa_sequence += 'X'

                ns1_aa_seq = Seq(aa_sequence)

                # Create SeqRecord for the amino acid sequence
                aa_record = SeqRecord(
                    ns1_aa_seq,
                    id=sample_name,
                    description=f"NS1 amino acid sequence from {sample_name}"
                )
                ns1_aa_records.append(aa_record)

                print(f"  NS1 amino acid length: {len(ns1_aa_seq)} aa")
                print(f"  First 20 aa: {str(ns1_aa_seq)[:20]}")

            except Exception as e:
                print(f"  Error translating sequence for {sample_name}: {e}")

        except Exception as e:
            print(f"  Error processing {fasta_file}: {e}")

    # Create output directory if it doesn't exist
    os.makedirs("phylogenetics", exist_ok=True)

    # Write all NS1 amino acid sequences to output file
    output_file = "phylogenetics/vsp_nas1_aa_sequences.fasta"
    print(f"\nWriting {len(ns1_aa_records)} sequences to {output_file}")

    with open(output_file, "w") as out_handle:
        SeqIO.write(ns1_aa_records, out_handle, "fasta")

    print(f"Successfully wrote NS1 amino acid sequences to {output_file}")

    # Print summary
    print(f"\nSummary:")
    print(f"  Processed {len(consensus_files)} consensus files")
    print(f"  Generated {len(ns1_aa_records)} NS1 amino acid sequences")
    print(f"  Output file: {output_file}")

if __name__ == "__main__":
    extract_ns1_aa_sequences()