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

## Benchmarking

This assembly workflow was benchmarked against a previous version of the assembly workflow which consisted of Trimmomatic + SPAdes v3.11.1. Organisms relevant for the NRLBM were selected from the [FDA-ARGOS BioProject](https://www.ncbi.nlm.nih.gov/bioproject/231221) (*Neisseria* spp., *Streptococcus* spp., *Staphylococcus* spp., and *Escherichia coli*). Draft assemblies produced using this workflow were compared with the reference genomes provided with the FDA-ARGOS project.

## Instructions (Dutch)

Deze pipeline is geïnstalleerden klaar om te gebruiken op het Lisa account van het NRLBM.

In het kort:

- Upload de data naar de map "input" via het WinSCP programma.
- Als de upload voltooid is, klik bovenin op "Opdrachten" en vervolgens op "Terminal openen" (kan ook d.m.v. Ctrl + Shift + T).
- Type "start" en druk op Enter. Dit moet zonder hoofdletters gespeld zijn.
- Er verschijnt tekst die begint met "Submitted batch job" gevolgd door een nummer. Dit betekent dat de analyse in de wacht is gezet. Er komt een mail binnen zodra de analyse voltooid is. De datum en tijd van het indienen wordt het volgnummer voor de analyse, in het format jaar.maand.dag_uur.minuut.seconde. Een analyse ingediend op 1 oktober 2021 om 16:24:19 krijgt dus het volgnummer `2021.10.01_16.24.19`.
- Na notificatie verschijnt er in de map "output" een map met het volgnummer. Deze map bevat de geassembleerde genomen (map "genomes"), een map met kwaliteitscontroles (map "qc_reports") en een Excel bestand "summary.xlsx".

Een uitgebreide handleiding is te vinden op de G-schijf van MMI-NRLBM.
