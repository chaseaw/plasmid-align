#!/usr/bin/env python

"""
Created November 2025
Authors: ChatGPT-4o and chasew

General pipeline for sequencing plasmids or linear amplicons using either 
Illumina paired-end data or Nanopore long reads.

Default behavior:
	• If a reference FASTA is provided: perform reference-based alignment
	• If no reference is provided: user must explicitly request de novo assembly

Supported inputs:
	• Illumina paired-end reads:   --illumina R1.fastq.gz R2.fastq.gz
	• Nanopore long reads:         --nanopore reads.fastq

Optional modes:
	• Reference alignment (default when --reference is given)
	• De novo assembly             (--de_novo)

Example usages:
	# Illumina reads aligned to a plasmid or linear amplicon reference
	python plasmid_seq.py --illumina sample_R1.fq.gz sample_R2.fq.gz \
						  --reference plasmid.fasta \
						  --output plasmid_align

	# Nanopore reads aligned to a reference (plasmid or amplicon)
	python plasmid_seq.py --nanopore nanopore.fastq \
						  --reference plasmid.fasta \
						  --output ont_align

	# Illumina de novo assembly
	python plasmid_seq.py --illumina R1.fq.gz R2.fq.gz \
						  --de_novo \
						  --output illumina_assembly

	# Nanopore de novo assembly
	python plasmid_seq.py --nanopore longreads.fq \
						  --de_novo \
						  --output ont_assembly
"""
import argparse
import subprocess
import os
import sys
import re
import shutil

def strip_fastq_extensions(filename):
	for ext in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
		if filename.endswith(ext):
			return filename.removesuffix(ext)
	return filename

def strip_from_R1_R2(name):
	# Strip from last occurrence of R1 or R2 (not case-sensitive)
	match = re.search(r"(.*?)([._-]?R[12])(?!.*R[12]).*", name, re.IGNORECASE)
	return match.group(1) if match else name

def run_command(cmd, description):
	print(f"\n[RUNNING] {description}: {' '.join(cmd)}")
	result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	# Don't print stdout unless explicitly needed
	if result.stderr:
		print(result.stderr.decode(), file=sys.stderr)

	if result.returncode != 0:
		print(f"[ERROR] {description} failed with return code {result.returncode}")
		sys.exit(result.returncode)

	return result

def run_fastp(reads, platform, output_dir, threads=2):
	print("\n[INFO] Running fastp for quality and adapter trimming...")
	R1 = os.path.basename(strip_fastq_extensions(reads[0]))
	file_prefix = strip_from_R1_R2(R1)
	if platform == "illumina":
		R2 = os.path.basename(strip_fastq_extensions(reads[1]))
		r1_trim = os.path.join(output_dir, f"{R1}.trimmed.fastq.gz")
		r2_trim = os.path.join(output_dir, f"{R2}.trimmed.fastq.gz")
		cmd = [
			"fastp", "-i", reads[0], "-I", reads[1],
			"-o", r1_trim, "-O", r2_trim,
			"-w", str(threads),
			"-h", os.path.join(output_dir, f"{file_prefix}_fastp.html"),
			"-j", os.path.join(output_dir, f"{file_prefix}_fastp.json"),
		]
		run_command(cmd, "fastp trimming (Illumina)")
		return [r1_trim, r2_trim]
	
	else:  # Nanopore
		trimmed = os.path.join(output_dir, f"{R1}.trimmed.fastq.gz")
		cmd = [
			"fastp", "-i", reads[0], "-o", trimmed,
			"-w", str(threads),
			"-h", os.path.join(output_dir, f"{file_prefix}_fastp.html"),
			"-j", os.path.join(output_dir, f"{file_prefix}_fastp.json"),
		]
		run_command(cmd, "fastp trimming (Nanopore)")
		return [trimmed]

def align_with_reference(reads, ref_fasta, platform, output_dir, threads=2):
	R1 = os.path.basename(strip_fastq_extensions(reads[0]))
	file_prefix = strip_from_R1_R2(R1)

	ref_copy = os.path.join(output_dir, os.path.basename(ref_fasta))
	shutil.copy2(ref_fasta, ref_copy)
	
	# Index the reference
	if platform == "illumina":
		run_command(["bwa", "index", ref_copy], "Index reference with BWA")
		cmd = ["bwa", "mem", "-t", str(threads), ref_copy, reads[0], reads[1]]
	else:
		run_command(["minimap2", "-d", f"{ref_copy}.mmi", ref_copy], "Index reference with minimap2")
		cmd = ["minimap2", "-t", str(threads), "-ax", "map-ont", f"{ref_copy}.mmi", reads[0]]

	sam_file = os.path.join(output_dir, f"{file_prefix}.sam")
	with open(sam_file, "w") as f:
		result = run_command(cmd, f"Align reads to reference ({platform})")
		f.write(result.stdout.decode())

	bam_file = os.path.join(output_dir, f"{file_prefix}.bam")

	run_command(["samtools", "view", "-@", str(threads), "-bS", sam_file, "-o", bam_file], "Convert SAM to BAM")
	os.remove(sam_file)
	run_command(["samtools", "sort", "-@", str(threads), bam_file, "-o", os.path.join(output_dir, f"{file_prefix}.sorted.bam")], "Sort BAM")
	os.remove(bam_file)
	run_command(["samtools", "index", os.path.join(output_dir, f"{file_prefix}.sorted.bam")], "Index BAM")
	print(f"\n[OUTPUT] Sorted BAM: {os.path.join(output_dir, f'{file_prefix}.sorted.bam')}")

def assemble_de_novo(reads, platform, output_dir, threads=2, genome_size="15k"):

	R1 = os.path.basename(strip_fastq_extensions(reads[0]))
	file_prefix = strip_from_R1_R2(R1)

	if platform == "illumina":
		cmd = ["spades.py", "-1", reads[0], "-2", reads[1],
			   "-t", str(threads), "-o", output_dir]
		run_command(cmd, "De novo assembly with SPAdes")
		final = os.path.join(output_dir, "contigs.fasta")

	elif platform == "nanopore":
		print("\n[INFO] Performing manual downsampling for Nanopore de novo assembly...")

		trimmed_fastq = reads[0]

		# Determine number of reads in the FASTQ
		print("[INFO] Counting reads in trimmed FASTQ...")
		count_cmd = ["gunzip", "-c", trimmed_fastq] if trimmed_fastq.endswith(".gz") else ["cat", trimmed_fastq]
		wc_output = subprocess.check_output(count_cmd)
		total_lines = wc_output.decode().count("\n")
		total_reads = total_lines // 4

		print(f"[INFO] Total reads available: {total_reads}")

		# Target ~200x coverage for small plasmids:
		# For ~15 kb and ~9.5 kb reads, ~300–400 reads = ~200x coverage
		target_reads = 400 if total_reads > 1000 else total_reads
		print(f"[INFO] Targeting {target_reads} reads (~200x coverage)")

		downsampled_fastq = os.path.join(output_dir, f"{file_prefix}.downsampled.fastq")

		# Probability fraction needed
		frac = min(1.0, target_reads / total_reads)
		print(f"[INFO] Using seqtk to sample fraction {frac:.4f}")

		cmd = ["seqtk", "sample", "-s42", trimmed_fastq, str(frac)]
		with open(downsampled_fastq, "w") as f:
			result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

		if result.returncode != 0:
			print(result.stderr.decode(), file=sys.stderr)
			sys.exit("[ERROR] seqtk downsampling failed.")

		print(f"[INFO] Downsampled FASTQ written to: {downsampled_fastq}")

		cmd = ["flye", "--nano-raw", downsampled_fastq,	"--out-dir", output_dir, "--genome-size", genome_size, "--threads", str(threads)]

		print("[INFO] Running Flye (manual downsampling mode)...")
		run_command(cmd, "De novo assembly with Flye")

		final = os.path.join(output_dir, "assembly.fasta")

	# Final output management
	final_named = os.path.join(output_dir, f"{file_prefix}.assembly.fa")
	shutil.copy2(final, final_named)
	os.remove(final)

	print(f"\n[OUTPUT] Assembled FASTA: {final_named}")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Plasmid sequencing analysis pipeline")
	parser.add_argument("--illumina", nargs=2, metavar=('R1', 'R2'), help="Illumina paired-end reads")
	parser.add_argument("--nanopore", metavar='ONT_READS', help="Nanopore fastq file")
	parser.add_argument("--reference", metavar='FASTA', help="Optional reference FASTA for alignment")
	parser.add_argument("--output", default="plasmid_seq_out", help="Output directory for all results (default: plasmid_seq_out)")
	parser.add_argument("--de_novo", action="store_true", help="Force de novo assembly (cannot be used with reference alignment)")
	parser.add_argument("--threads", type=int, default=2, help="Number of threads to use (default: 2)")
	parser.add_argument("--genome_size", default="15k", help="Expected genome size for Flye (default: 15k for plasmids)")
	args = parser.parse_args()

	if args.illumina:
		platform = "illumina"
		reads = args.illumina
		assemble_suffix = "spades"
	elif args.nanopore:
		platform = "nanopore"
		reads = [args.nanopore]
		assemble_suffix = "flye"
	else:
		print("[ERROR] Please provide either --illumina or --nanopore reads.")
		sys.exit(1)

	if args.illumina and args.nanopore:
		print("[ERROR] Please specify either --illumina or --nanopore, not both.")
		sys.exit(1)

	if args.reference and args.de_novo:
		print("[ERROR] --reference and --de_novo cannot be used together.")
		sys.exit(1)

	output_dir = os.path.abspath(args.output)
	trim_dir = os.path.join(output_dir, "trimmed")
	align_or_assemble_dir = os.path.join(output_dir, "alignment" if args.reference else f"assembly_{assemble_suffix}")

	os.makedirs(output_dir, exist_ok=True)
	os.makedirs(trim_dir, exist_ok=True)
	os.makedirs(align_or_assemble_dir, exist_ok=True)
	sys.stdout = open(os.path.join(output_dir, "pipeline.out"), "w")
	sys.stderr = sys.stdout

	reads = run_fastp(reads, platform, trim_dir, args.threads)

	if args.reference:
		align_with_reference(reads, args.reference, platform, align_or_assemble_dir, threads=args.threads)
	elif args.de_novo:
		assemble_de_novo(reads, platform, align_or_assemble_dir, threads=args.threads, genome_size=args.genome_size)
	else:
		print("[ERROR] No reference provided and --de_novo flag not set. Exiting.")
		sys.exit(1)