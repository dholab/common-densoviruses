#!/usr/bin/env python3
"""
Script to prepare SRA paired-end input files for alignment
Downloads SRA and keeps R1/R2 separate (no interleaving or compression)
This is faster and more memory efficient than concatenating/compressing
"""

import os
import subprocess


def get_num_threads():
    """Get number of threads from snakemake or default to 4"""
    try:
        return snakemake.threads
    except:
        return 4


def download_sra_paired(sra_accession, output_r1, output_r2, threads=4):
    """Download SRA accession and keep R1/R2 files separate"""
    print(f"Downloading SRA accession (paired-end): {sra_accession}")

    # Get output directory
    output_dir = os.path.dirname(output_r1)
    if not output_dir:
        output_dir = '.'

    os.makedirs(output_dir, exist_ok=True)

    # Use fasterq-dump.3.1.1 (available in the apptainer)
    cmd = [
        'fasterq-dump.3.1.1',
        sra_accession,
        '-O', output_dir,
        '-e', '1'  # Use 1 thread for download
    ]

    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # fasterq-dump outputs uncompressed files as {accession}_1.fastq and {accession}_2.fastq
    fastq_1 = os.path.join(output_dir, f"{sra_accession}_1.fastq")
    fastq_2 = os.path.join(output_dir, f"{sra_accession}_2.fastq")

    # Check that both files exist
    if not os.path.exists(fastq_1) or not os.path.exists(fastq_2):
        raise FileNotFoundError(
            f"Could not find paired-end output FASTQ files from fasterq-dump. "
            f"Expected {fastq_1} and {fastq_2}"
        )

    # Rename to match expected output filenames (based on sample name, not SRA accession)
    if fastq_1 != output_r1:
        os.rename(fastq_1, output_r1)
    if fastq_2 != output_r2:
        os.rename(fastq_2, output_r2)

    print(f"Successfully prepared paired-end files:")
    print(f"  R1: {output_r1}")
    print(f"  R2: {output_r2}")


def main():
    # Get parameters from snakemake
    output_r1 = snakemake.output.r1
    output_r2 = snakemake.output.r2
    input_type = snakemake.params.input_type
    input_files = snakemake.params.input_files
    threads = get_num_threads()

    print(f"Preparing SRA paired-end input")
    print(f"Input type: {input_type}")

    # This script only handles SRA paired-end
    if input_type != 'sra_paired':
        raise ValueError(f"This script only handles sra_paired input type, got: {input_type}")

    download_sra_paired(input_files['sra'], output_r1, output_r2, threads)


if __name__ == '__main__':
    main()
