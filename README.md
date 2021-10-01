# Bacterial isolate genome assembly workflow for NRLBM

## Introduction

This repository contains code for the basic assembly workflow used by the NRLBM. The workflow contains assembly using SPAdes and some basic quality controls around it. The pipeline is written in Snakemake and is submitted to a SLURM scheduler.

## Methods

The pipeline consists of several steps:

- `workflow/scripts/check.py`: Python script that checks a couple of things before the job is submitted to the SLURM scheduler. This is the only step not submitted to the cluster.
- FastQC (before and after trimming): checks the read set quality
- Trimmomatic: trims and filters reads 
- Shovill: assembles trimmed reads
- Quast: evaluates assembly quality
- mlst: assigns an MLST profile based on the assembled genome
- Kraken2: identifies species present in trimmed reads
- MultiQC: collects quality control reports for comparison
- coverage: calculates sequencing read depth
