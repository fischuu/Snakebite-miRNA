# Snakebite-miRNA

A Snakemake pipeline for small RNA-seq / miRNA-seq analysis: adapter trimming, tRNA/PhiX
decontamination, alignment against species-filtered miRBase mature/hairpin references and
the host reference genome (Bowtie and STAR), miRNA and gene-level quantification, and
novel-miRNA locus prediction from genome coverage. It runs on HPC (SLURM) with every tool
delivered as a Singularity/Apptainer container, and follows the same module/escalation/
report architecture as [Snakebite-Holoruminant-MetaG](https://github.com/fischuu/Snakebite-Holoruminant-MetaG).

# Installation

```
git clone git@github.com:fischuu/Snakebite-miRNA.git
```

Or via HTTPS:

```
git clone https://github.com/fischuu/Snakebite-miRNA.git
```

# Running the pipeline

Copy `run_Snakebite-miRNA.sh` into your project folder, point `projectFolder` at your
project and `pipeline_folder` (in `config/config.yaml`) at this checkout, then:

```bash
bash run_Snakebite-miRNA.sh            # full pipeline
bash run_Snakebite-miRNA.sh align -np  # dry run of a single module
```

Modules (top-level rules in `rule all`): `reads`, `reference`, `preprocess`,
`decontaminate`, `align`, `novel_mirna`, `quantify`. Each module can be run as a whole or
as a single tool within it using `<module>__<tool>` naming, and a `report_<module>` target
renders that module's PDF report.

# Included tools

 * Bowtie
 * STAR
 * cutadapt
 * FastQC / MultiQC
 * samtools
 * subread (featureCounts)
 * seqkit
 * bedtools
 * extractSoftclipped (vendored)

# Configuration

 * `config/config.yaml` -- paths to the other config files and the SLURM resource tiers
   (`resource_sets`).
 * `config/features.yaml` -- paths to the tRNA/PhiX/mature/hairpin/genome references and
   the genome annotation.
 * `config/params.yaml` -- species identifier, adapter-trimming protocol/preset, and every
   other tool parameter.
 * `config/escalation.yaml` -- per-rule ordered list of resource tiers to escalate through
   on retry.
 * `config/samples.tsv` -- sample sheet: `sample_id`, `library_id` (one row per sequencing
   lane), `read1`, and the `group`/`de_group` columns used for downstream differential
   analysis in the report.
 * `config/.docker.yml` -- one container image per tool.
 * `config/profiles/<cluster>/config.yaml` -- SLURM executor profile; copy `Puhti/` for a
   new cluster.
