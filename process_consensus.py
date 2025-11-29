#!/usr/bin/env python3

import sys
import os
import re

def get_logan_metadata(metadata_file):
    metadata = {}
    with open(metadata_file, 'r') as f:
        # Skip header
        for _ in range(4):
            next(f)
        for line in f:
            if not line.strip() or line.startswith('|--'):
                continue
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 7:
                biosample = parts[1]
                organism = parts[6]
                # Extract the biosample accession number
                match = re.search(r'\[(SAMN[0-9]+)\]', biosample)
                if match:
                    accession = match.group(1)
                    metadata[accession] = organism.replace(' ', '_')
    return metadata

def translate_dna(dna):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if len(codon) == 3:
            protein += codon_table.get(codon.upper(), 'X')
    return protein

def process_fasta(file_path, metadata):
    with open(file_path, 'r') as f:
        header = ''
        sequence = ''
        for line in f:
            if line.startswith('>'):
                header = line.strip()
            else:
                sequence += line.strip()

    if not sequence:
        return

    ns1_part1 = sequence[237:301]
    ns1_part2 = sequence[372:2015]
    ns1_dna = ns1_part1 + ns1_part2

    if len(ns1_dna) != 1707:
        # sys.stderr.write(f"Warning: NS1 sequence in {file_path} is not complete. Length is {len(ns1_dna)}. Skipping.\n")
        return

    ns1_protein = translate_dna(ns1_dna)

    base_name = os.path.basename(file_path)
    sample_name = base_name.replace('.consensus.fasta', '')

    organism = "Unknown"
    if "logan/consensus" in file_path:
        organism = metadata.get(sample_name, "Unknown_species")
    elif "air-samples/vsp/consensus" in file_path:
        organism = "Human_associated"


    print(f">{sample_name}_{organism}_NS1")
    print(ns1_protein)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: ./process_consensus.py <metadata_file> <fasta_files...>")
        sys.exit(1)
    
    metadata_file = sys.argv[1]
    fasta_files = sys.argv[2:]

    metadata = get_logan_metadata(metadata_file)

    for fasta_file in fasta_files:
        process_fasta(fasta_file, metadata)
