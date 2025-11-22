#!/usr/bin/env python3
"""
Script to prepare input files for alignment
Handles SRA downloads (single-end or interleaved) and single-end FASTQ
Supports multiple SRA accessions per sample (semicolon-separated)
Uses pigz for parallel compression (much faster than Python's gzip)

Note: For SRA paired-end reads that should be kept separate, use prepare_input_paired.py
"""

import os
import sys
import subprocess
import shutil

def get_num_threads():
    """Get number of threads from snakemake or default to 4"""
    try:
        return snakemake.threads
    except:
        return 4


def download_single_sra(sra_accession, output_dir, threads=4):
    """Download a single SRA accession and return the paths to the output files.

    Returns:
        tuple: (r1_file, r2_file) for paired-end, or (single_file, None) for single-end
    """
    print(f"Downloading SRA accession: {sra_accession}")

    # Use fasterq-dump.3.1.1 (available in the apptainer)
    cmd = [
        'fasterq-dump.3.1.1',
        sra_accession,
        '-O', output_dir,
        '-e', '1'  # Use 1 thread for download
    ]

    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # fasterq-dump outputs uncompressed files
    # For paired-end, it creates _1.fastq and _2.fastq
    # For single-end, it creates .fastq

    fastq_1 = os.path.join(output_dir, f"{sra_accession}_1.fastq")
    fastq_2 = os.path.join(output_dir, f"{sra_accession}_2.fastq")
    fastq_single = os.path.join(output_dir, f"{sra_accession}.fastq")

    if os.path.exists(fastq_1) and os.path.exists(fastq_2):
        return (fastq_1, fastq_2)
    elif os.path.exists(fastq_single):
        return (fastq_single, None)
    else:
        raise FileNotFoundError(
            f"Could not find output FASTQ files from fasterq-dump. "
            f"Expected {fastq_1}/{fastq_2} or {fastq_single}"
        )


def download_multiple_sra(sra_accessions, output_fastq, threads=4):
    """Download multiple SRA accessions and concatenate them.

    For paired-end data: concatenates all R1 files and all R2 files separately
    For single-end data: concatenates all files together

    Args:
        sra_accessions: List of SRA accession numbers
        output_fastq: Path to output fastq.gz file (or marker file for paired-end)
        threads: Number of threads for compression
    """
    print(f"Downloading {len(sra_accessions)} SRA accessions: {', '.join(sra_accessions)}")

    output_dir = os.path.dirname(output_fastq)
    if not output_dir:
        output_dir = '.'

    os.makedirs(output_dir, exist_ok=True)

    # Download all accessions
    all_r1_files = []
    all_r2_files = []
    all_single_files = []

    for sra_acc in sra_accessions:
        r1, r2 = download_single_sra(sra_acc, output_dir, threads)
        if r2 is not None:
            all_r1_files.append(r1)
            all_r2_files.append(r2)
        else:
            all_single_files.append(r1)

    # Get sample name from output path
    sample_name = os.path.basename(output_fastq).replace('.fastq.gz', '')

    # Handle the downloaded files
    if all_r1_files and all_r2_files:
        # Paired-end data - concatenate R1 and R2 files separately
        print(f"Detected paired-end reads from {len(all_r1_files)} accessions, concatenating...")

        output_r1 = os.path.join(output_dir, f"{sample_name}_1.fastq")
        output_r2 = os.path.join(output_dir, f"{sample_name}_2.fastq")

        # Concatenate R1 files
        if len(all_r1_files) == 1:
            os.rename(all_r1_files[0], output_r1)
        else:
            r1_files_str = " ".join(all_r1_files)
            cmd = f"cat {r1_files_str} > {output_r1}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            # Clean up individual files
            for f in all_r1_files:
                os.remove(f)

        # Concatenate R2 files
        if len(all_r2_files) == 1:
            os.rename(all_r2_files[0], output_r2)
        else:
            r2_files_str = " ".join(all_r2_files)
            cmd = f"cat {r2_files_str} > {output_r2}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            # Clean up individual files
            for f in all_r2_files:
                os.remove(f)

        # Create marker file for snakemake
        with open(output_fastq, 'w') as f:
            f.write(f"# Paired-end data from {len(sra_accessions)} SRA accessions\n")
            f.write(f"# Accessions: {', '.join(sra_accessions)}\n")
            f.write(f"# See {sample_name}_1.fastq and {sample_name}_2.fastq\n")

        print(f"Successfully prepared paired-end files from {len(sra_accessions)} accessions:")
        print(f"  R1: {output_r1}")
        print(f"  R2: {output_r2}")

    elif all_single_files:
        # Single-end data - concatenate and compress
        print(f"Detected single-end reads from {len(all_single_files)} accessions, concatenating and compressing...")

        if len(all_single_files) == 1:
            # Just compress the single file
            cmd = f"pigz -p {threads} -c {all_single_files[0]} > {output_fastq}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            os.remove(all_single_files[0])
        else:
            # Concatenate and compress
            single_files_str = " ".join(all_single_files)
            cmd = f"cat {single_files_str} | pigz -p {threads} > {output_fastq}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            # Clean up individual files
            for f in all_single_files:
                os.remove(f)

        print(f"Successfully prepared: {output_fastq}")

    else:
        raise ValueError("No FASTQ files were downloaded successfully")


def download_sra(sra_input, output_fastq, threads=4):
    """Download SRA accession(s) and convert to FASTQ.

    Args:
        sra_input: Single SRA accession or semicolon-separated list of accessions
        output_fastq: Path to output file
        threads: Number of threads for compression
    """
    # Handle multiple SRA accessions (semicolon-separated)
    if ';' in sra_input:
        sra_accessions = [s.strip() for s in sra_input.split(';') if s.strip()]
        download_multiple_sra(sra_accessions, output_fastq, threads)
    else:
        # Single accession - use the multi-SRA function for consistency
        download_multiple_sra([sra_input.strip()], output_fastq, threads)


def get_sra_accessions_for_biosample(biosample_id):
    """Fetch all SRA run accessions associated with a BioSample ID.

    Uses NCBI E-utilities (esearch/efetch) to query the SRA database.

    Args:
        biosample_id: BioSample accession (e.g., SAMN12345678)

    Returns:
        List of SRA run accessions (e.g., ['SRR1234567', 'SRR1234568'])
    """
    print(f"Fetching SRA accessions for BioSample: {biosample_id}")

    # Use esearch and efetch to get run info
    # esearch -db sra -query "SAMN12345678" | efetch -format runinfo
    cmd = f'esearch -db sra -query "{biosample_id}" | efetch -format runinfo'

    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error fetching SRA accessions: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"Failed to fetch SRA accessions for {biosample_id}")

    # Parse the runinfo CSV output - first column is Run accession
    # Skip header line and extract run accessions
    sra_accessions = []
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        # First column is the Run accession
        fields = line.split(',')
        if fields and fields[0]:
            acc = fields[0].strip()
            # Filter for valid SRA/ERA/DRA accessions
            if acc.startswith(('SRR', 'ERR', 'DRR')):
                sra_accessions.append(acc)

    # Remove duplicates and sort
    sra_accessions = sorted(set(sra_accessions))

    if not sra_accessions:
        raise ValueError(f"No SRA accessions found for BioSample {biosample_id}")

    print(f"Found {len(sra_accessions)} SRA accessions: {', '.join(sra_accessions)}")
    return sra_accessions


def download_biosample(biosample_id, output_fastq, threads=4):
    """Download all SRA accessions for a BioSample and prepare FASTQ files.

    Args:
        biosample_id: BioSample accession (e.g., SAMN12345678)
        output_fastq: Path to output file
        threads: Number of threads for compression
    """
    # Get all SRA accessions for this biosample
    sra_accessions = get_sra_accessions_for_biosample(biosample_id)

    # Download all accessions
    download_multiple_sra(sra_accessions, output_fastq, threads)

def prepare_paired_fastq(r1_file, r2_file, output_fastq, threads=4):
    """Concatenate paired-end FASTQ files using pigz for fast compression"""
    print(f"Concatenating paired-end files: {r1_file} and {r2_file}")

    # Build command based on whether inputs are gzipped
    if r1_file.endswith('.gz') and r2_file.endswith('.gz'):
        # Both gzipped: decompress with pigz, concatenate, recompress
        # Note: This is still slower than passing directly to minimap2,
        # but much faster than Python's gzip
        cmd = f"(pigz -dc -p {threads} {r1_file}; pigz -dc -p {threads} {r2_file}) | pigz -p {threads} > {output_fastq}"
    elif r1_file.endswith('.gz'):
        cmd = f"(pigz -dc -p {threads} {r1_file}; cat {r2_file}) | pigz -p {threads} > {output_fastq}"
    elif r2_file.endswith('.gz'):
        cmd = f"(cat {r1_file}; pigz -dc -p {threads} {r2_file}) | pigz -p {threads} > {output_fastq}"
    else:
        cmd = f"cat {r1_file} {r2_file} | pigz -p {threads} > {output_fastq}"

    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"Successfully prepared: {output_fastq}")

def prepare_single_fastq(input_file, output_fastq, threads=4):
    """Prepare single-end FASTQ file"""
    print(f"Preparing single-end file: {input_file}")

    if input_file.endswith('.gz'):
        # Already compressed, just copy (or symlink for speed)
        shutil.copy2(input_file, output_fastq)
    else:
        # Compress with pigz
        cmd = f"pigz -p {threads} -c {input_file} > {output_fastq}"
        print(f"Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    print(f"Successfully prepared: {output_fastq}")

def main():
    # Get parameters from snakemake
    output_fastq = snakemake.output.fastq
    input_type = snakemake.params.input_type
    input_files = snakemake.params.input_files
    threads = get_num_threads()

    print(f"Using {threads} threads for compression")

    # Create output directory if needed
    os.makedirs(os.path.dirname(output_fastq), exist_ok=True)

    # Process based on input type
    if input_type == 'sra':
        download_sra(input_files['sra'], output_fastq, threads)
    elif input_type == 'biosample':
        download_biosample(input_files['biosample'], output_fastq, threads)
    elif input_type == 'paired':
        prepare_paired_fastq(input_files['r1'], input_files['r2'], output_fastq, threads)
    elif input_type == 'single':
        prepare_single_fastq(input_files['se'], output_fastq, threads)
    else:
        raise ValueError(f"Unknown input type: {input_type}")

if __name__ == '__main__':
    main()
