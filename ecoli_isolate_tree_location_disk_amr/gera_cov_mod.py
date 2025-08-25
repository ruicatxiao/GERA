#!/usr/bin/env python3
import os
import pandas as pd
import subprocess
import logging
import argparse
import sys
from pathlib import Path

def setup_logging(output_dir):
    """Configure logging to file and console"""
    log_path = output_dir / "logs" / "pipeline.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # File handler
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger

def run_command(cmd, sample, logger):
    """Execute system command with verbose output handling"""
    try:
        process = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )

        # Stream output in real-time
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                logger.info(f"[{sample}] {output.strip()}")

        exit_code = process.poll()
        if exit_code != 0:
            raise subprocess.CalledProcessError(exit_code, cmd)
            
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed for {sample} (exit {e.returncode}): {e.cmd}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error for {sample}: {str(e)}")
        return False

def process_sample(sample, barcode, bam_folder, args, fasta_path, logger):
    """Process a single sample through the pipeline"""
    try:
        # Create output directories
        (args.output/"bam").mkdir(exist_ok=True)
        (args.output/"mod").mkdir(exist_ok=True)

        # File paths
        bam_path = Path(bam_folder) / f"SQK-NBD114-24_{barcode}.bam"
        outputs = {
            'raw_bam': args.output/"bam"/f"{sample}.bam",
            'sorted_bam': args.output/"bam"/f"{sample}_sorted.bam",
            'coverage': args.output/"bam"/f"{sample}.txt",
            'bedmethyl': args.output/"mod"/f"{sample}.bedmethyl"
        }

        logger.info(f"Processing {sample} with BAM: {bam_path}")

        # 1. Dorado alignment with verbose output
        cmd = (
            f"dorado aligner --verbose --threads {args.threads} "
            f"{fasta_path} {bam_path} > {outputs['raw_bam']}"
        )
        if not run_command(cmd, sample, logger):
            return False

        # 2. Sort BAM
        cmd = f"samtools sort -@ {args.threads} -o {outputs['sorted_bam']} {outputs['raw_bam']}"
        if not run_command(cmd, sample, logger):
            return False

        # 3. Index BAM
        cmd = f"samtools index -@ {args.threads} {outputs['sorted_bam']}"
        if not run_command(cmd, sample, logger):
            return False

        # 4. Coverage analysis
        cmd = f"samtools coverage -d 0 {outputs['sorted_bam']} > {outputs['coverage']}"
        if not run_command(cmd, sample, logger):
            return False

        # 5. Modkit pileup
        cmd = f"modkit pileup -t {args.threads} --filter-threshold 0.95 {outputs['sorted_bam']} {outputs['bedmethyl']}"
        if not run_command(cmd, sample, logger):
            return False

        return True

    except Exception as e:
        logger.error(f"Critical error processing {sample}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Simplified Methylation Pipeline")
    parser.add_argument("--fa_input", required=True, help="Directory with FASTA files")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--sample_barcode", required=True, help="Path to sample-barcode lookup table")
    parser.add_argument("--threads", type=int, default=32)
    args = parser.parse_args()

    # Validate inputs
    if not Path(args.fa_input).exists():
        sys.exit(f"Error: FASTA directory not found: {args.fa_input}")
    args.output.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger = setup_logging(args.output)
    
    # Load sample data
    try:
        samples_df = pd.read_csv(args.sample_barcode, sep='\t')
        logger.info(f"Loaded {len(samples_df)} samples from lookup table")
    except Exception as e:
        sys.exit(f"Error loading sample table: {str(e)}")

    # Map FASTA files
    fasta_files = {}
    for fasta in Path(args.fa_input).glob("*.fasta"):
        sample_id = fasta.name.split('_')[0]
        fasta_files[sample_id] = fasta.resolve()
    logger.info(f"Mapped {len(fasta_files)} FASTA files")

    # Process samples
    success_count = 0
    for _, row in samples_df.iterrows():
        sample = row['sample']
        barcode = row['barcode']
        bam_folder = row['folder']

        if sample not in fasta_files:
            logger.error(f"Skipping {sample} - no FASTA file")
            continue

        success = process_sample(
            sample=sample,
            barcode=barcode,
            bam_folder=bam_folder,
            args=args,
            fasta_path=fasta_files[sample],
            logger=logger
        )

        if success:
            success_count += 1
            logger.info(f"Successfully processed {sample}")
        else:
            logger.error(f"Failed to process {sample}")

    logger.info(f"Pipeline complete. Successfully processed {success_count}/{len(samples_df)} samples")

if __name__ == "__main__":
    main()