#!/usr/bin/env python3
"""
Extract NS1 amino acid sequences from Logan consensus files.
Uses metadata from logan/sra-metadata.md to create proper FASTA headers.
"""

import os
import re
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_metadata():
    """Parse the SRA metadata markdown table to create BioSample -> Organism mapping."""
    metadata_file = "logan/sra-metadata.md"
    biosample_to_organism = {}

    print("Parsing metadata from", metadata_file)

    with open(metadata_file, 'r') as f:
        lines = f.readlines()

    # Find the header line and data lines
    header_found = False
    for line in lines:
        line = line.strip()

        # Skip until we find the table header
        if '| BioSample |' in line:
            header_found = True
            continue

        # Skip the separator line
        if header_found and line.startswith('|---'):
            continue

        # Process data lines
        if header_found and line.startswith('|') and '|' in line:
            # Split the line by pipes and clean up
            parts = [part.strip() for part in line.split('|')]

            if len(parts) >= 7:  # Make sure we have enough columns
                biosample_raw = parts[1]  # BioSample column (with markdown link)
                organism = parts[6]      # Organism column

                # Extract BioSample ID from markdown link format [SAMN...](url)
                # or just plain text
                biosample_match = re.search(r'\[([^\]]+)\]', biosample_raw)
                if biosample_match:
                    biosample = biosample_match.group(1)
                else:
                    biosample = biosample_raw.strip()

                # Skip if this is an error row
                if organism.upper() == 'ERROR' or biosample.upper() == 'ERROR':
                    continue

                # Clean up organism name
                organism = organism.strip()
                if organism and biosample:
                    biosample_to_organism[biosample] = organism

    print(f"Parsed metadata for {len(biosample_to_organism)} samples")
    return biosample_to_organism

def extract_logan_ns1_aa_sequences():
    """Extract NS1 amino acid sequences from all Logan consensus files."""

    # Parse metadata for naming
    biosample_to_organism = parse_metadata()

    # NS1 CDS coordinates from GFF file (1-based, converting to 0-based for Python)
    # CDS 238-301 and CDS 373-2015
    ns1_cds1_start = 238 - 1  # 237 (0-based)
    ns1_cds1_end = 301       # 301 (0-based, exclusive)
    ns1_cds2_start = 373 - 1  # 372 (0-based)
    ns1_cds2_end = 2015      # 2015 (0-based, exclusive)

    # Find all Logan consensus files
    consensus_files = glob.glob("logan/consensus/*.consensus.fasta")
    consensus_files.sort()

    # List to store all NS1 amino acid sequences
    ns1_aa_records = []

    print(f"Found {len(consensus_files)} Logan consensus files to process")

    for fasta_file in consensus_files:
        print(f"Processing {fasta_file}")

        # Extract BioSample from filename (e.g., SAMN03282914.consensus.fasta -> SAMN03282914)
        filename = os.path.basename(fasta_file)
        biosample = filename.replace(".consensus.fasta", "")

        # Get organism name from metadata
        organism = biosample_to_organism.get(biosample, "unknown_organism")

        # Create FASTA header: [BioSample]_[Organism] with spaces as underscores
        fasta_header = f"{biosample}_{organism.replace(' ', '_')}"

        print(f"  BioSample: {biosample}")
        print(f"  Organism: {organism}")
        print(f"  FASTA header: {fasta_header}")

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
                print(f"  Warning: NS1 CDS length ({len(ns1_complete_cds)}) not divisible by 3 for {biosample}")
                # Trim to the largest multiple of 3
                ns1_complete_cds = ns1_complete_cds[:len(ns1_complete_cds) - (len(ns1_complete_cds) % 3)]
                print(f"  Trimmed to {len(ns1_complete_cds)} bp")

            # Translate to amino acids, handling ambiguous nucleotides properly
            try:
                # Handle ambiguous nucleotides by translating codon by codon
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
                # Replace all whitespaces with underscores in the description as well
                description = f"NS1_amino_acid_sequence_from_{biosample}_({organism.replace(' ', '_')})"

                aa_record = SeqRecord(
                    ns1_aa_seq,
                    id=fasta_header,
                    description=description
                )
                ns1_aa_records.append(aa_record)

                print(f"  NS1 amino acid length: {len(ns1_aa_seq)} aa")
                print(f"  First 20 aa: {str(ns1_aa_seq)[:20]}")

            except Exception as e:
                print(f"  Error translating sequence for {biosample}: {e}")

        except Exception as e:
            print(f"  Error processing {fasta_file}: {e}")

    # Create output directory if it doesn't exist
    os.makedirs("phylogenetics", exist_ok=True)

    # Write all NS1 amino acid sequences to output file
    output_file = "phylogenetics/logan_ns1_aa_sequences.fasta"
    print(f"\nWriting {len(ns1_aa_records)} sequences to {output_file}")

    with open(output_file, "w") as out_handle:
        SeqIO.write(ns1_aa_records, out_handle, "fasta")

    print(f"Successfully wrote NS1 amino acid sequences to {output_file}")

    # Print summary
    print(f"\nSummary:")
    print(f"  Processed {len(consensus_files)} consensus files")
    print(f"  Generated {len(ns1_aa_records)} NS1 amino acid sequences")
    print(f"  Output file: {output_file}")

    # Show a few example headers
    print(f"\nExample FASTA headers:")
    for i, record in enumerate(ns1_aa_records[:5]):
        print(f"  {record.id}")
    if len(ns1_aa_records) > 5:
        print(f"  ... and {len(ns1_aa_records) - 5} more")

if __name__ == "__main__":
    extract_logan_ns1_aa_sequences()