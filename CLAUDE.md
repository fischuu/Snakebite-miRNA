# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Snakemake pipeline for small RNA-seq / miRNA-seq analysis (Luke, Natural Resources
Institute Finland): adapter trimming, tRNA/PhiX decontamination, alignment against
species-filtered miRBase mature/hairpin references and the host reference genome (Bowtie
and STAR), miRNA and gene-level quantification, and novel-miRNA locus prediction from
per-sample genome coverage. It follows the same module/escalation/docker/report
architecture as `Snakebite-Holoruminant-MetaG` (a sibling pipeline in the same "Snakebite"
family) — that repo's `CLAUDE.md` is a useful cross-reference for the shared conventions.
It runs on HPC (SLURM) with every tool delivered as a Singularity/Apptainer container.

## Running the pipeline

There is no build step and no test suite. All "development" happens by editing `.smk` rule
files and dry-running Snakemake.

The pipeline is invoked through the wrapper script, never by calling `snakemake` directly:

```bash
bash run_Snakebite-miRNA.sh <module_or_rule> [snakemake args...]

# dry run (always do this first when editing rules)
bash run_Snakebite-miRNA.sh <module_or_rule> -np
```

The wrapper hardcodes `projectFolder`/`configFile` paths at the top and calls `snakemake`
with `--use-singularity`, a SLURM `--profile`, and `--restart-times 0 --retries 0` (retry/
escalation is instead handled per-rule, see "Resource escalation" below). This is a
*project* wrapper — copy it into each deployment and point it at that project's paths.

Modules (top-level rules in `rule all`): `reads`, `reference`, `preprocess`,
`decontaminate`, `align`, `novel_mirna`, `quantify`. Each module can be run as a whole or as
a single tool within it using `<module>__<tool>` naming, e.g.:

```bash
bash run_Snakebite-miRNA.sh align__star__mature
bash run_Snakebite-miRNA.sh novel_mirna
```

Running a submodule automatically triggers its upstream dependencies. A `report_<module>`
target exists for every module and renders that module's PDF report.

## Architecture

### Config-file driven, not code-driven

Everything a user needs to change lives in `config/`, loaded once at the top of
`workflow/Snakefile`:

- `config/config.yaml` — paths to the other config files, `pipeline_folder`, temp-storage
  fallback chain, and `resource_sets` (named SLURM resource profiles: `small`, `medium`,
  `medium_nvme`, `large`, `large_nvme`, `highmem`, `longrun`, `highmem_longrun`,
  `full_node` — shared tier names with Snakebite-Holoruminant-MetaG, though this pipeline
  doesn't need every tier yet).
- `config/features.yaml` — paths to the tRNA/PhiX/mature/hairpin/genome references and the
  genome annotation.
- `config/params.yaml` — species identifier (`params["species_id"]`, used to filter
  mature.fa/hairpin.fa down to one species), the adapter-trimming protocol preset
  (`illumina`/`nextflex`/`qiagen`/`lexogen`, applied in `workflow/Snakefile` by overwriting
  `params["preprocess"]["cutadapt"]` — leave `protocol` empty to configure cutadapt
  yourself), and every other tool parameter.
- `config/escalation.yaml` — per-rule ordered list of `resource_sets` to escalate through
  on retry, keyed by full rule name. Rules not listed default to
  `["small", "medium", "large"]` (`get_escalation_order` in
  `workflow/rules/helpers/__functions__.smk`).
- `config/samples.tsv` — sample sheet: `sample_id`, `library_id` (one row per sequencing
  lane), `read1`, and `group`/`de_group` (merges what used to be two separate files,
  `sampleSheet.tsv` and `sampleInfo.tsv`, in the pre-rewrite pipeline).
- `config/.docker.yml` — one container image per tool (not per module — this pipeline mixes
  several distinct tools per module, e.g. `align` uses both `bowtie` and `star` and
  `samtools`).
- `config/profiles/<cluster>/config.yaml` — SLURM executor profile; copy `Puhti/` for a new
  cluster.

`workflow/Snakefile` parses `samples.tsv` into `SAMPLES` (unique `sample_id`s) and
`SAMPLE_LIBRARY` (`(sample_id, library_id)` pairs) globals used throughout the rule files.

### Module layout under `workflow/rules/`

Each module is a directory with `__main__.smk` (includes the module's rule files and
defines the aggregate `rule <module>:`), `__functions__.smk` (module-local wildcard/lookup
helpers), and one `.smk` file per tool. `workflow/Snakefile` includes modules in pipeline
order: `helpers` → `reads` → `reference` → `preprocess` → `decontaminate` → `align` →
`novel_mirna` → `quantify` → `report`. `novel_mirna` is included *before* `quantify` even
though it conceptually comes after alignment, because `quantify__bedtools__novel_mirna`
needs `novel_mirna`'s predicted-loci BED file as an input.

`workflow/rules/folders.smk` centralizes every output-path constant (`MATURE_BAM_STAR`,
`TRNA_FASTQ`, `MPILEUP`, etc., all `pathlib.Path` objects) — check here first when tracing
where a rule's output goes.

### The read-processing chain

Reads flow through the modules in a specific, mostly-linear chain (see
`workflow/rules/align/star.smk` and `bowtie.smk` for the exact wiring):

```
raw lane (read1) -> preprocess__cutadapt (per lane)
                  -> preprocess__concatenate (lanes -> one FASTQ per sample)
                  -> decontaminate__bowtie__trna    (tRNA-unmapped reads continue)
                  -> decontaminate__bowtie__phix    (PhiX-unmapped reads continue)
                  -> align__bowtie__mature -> align__bowtie__hairpin -> align__bowtie__genome (kept, not a default target)
                  -> align__star__mature        (parallel branch off PhiX-unmapped)
                  -> align__star__mature_species -> align__star__hairpin_species -> align__star__genome
                  -> align__star__hairpin       (parallel branch off PhiX-unmapped, feeds mirbase quantification only)
```

The STAR `mature`/`hairpin` branches and the `mature_species` -> `hairpin_species` ->
`genome` chain are **both** rooted at PhiX-unmapped reads but are otherwise independent —
don't assume `align__star__hairpin`'s output feeds `align__star__genome`; it doesn't. Only
the `_species` branch does.

### Resource escalation (retry-with-bigger-resources)

Identical mechanism to Snakebite-Holoruminant-MetaG: every rule uses
`esc(key, rule_name)` / `esc_val(key, rule_name)` (from
`workflow/rules/helpers/__functions__.smk`) instead of hardcoding `resources:`/`threads:`.
`retries:` on each rule is `len(get_escalation_order(rule_name))`, and a loop at the bottom
of `workflow/Snakefile` forces `r.retries = 0` for any rule whose escalation list has ≤1
entry. Follow the existing pattern when adding a new rule, and add an entry to
`escalation.yaml` if it needs more than the default 3-tier escalation.

### Report module

`workflow/scripts/R/report_helpers.R` holds shared runtime/resource-allocation helpers and
tool-output parsers, sourced by every per-module `.Rmd` in the same folder. Benchmark files
live under a flat `benchmark/<module>/<rule_slug>.<wildcards...>.tsv` convention at the
project root (not nested under a `results/` tree) — `discover_benchmarks()` derives a
rule's report label as everything before the first `.` in the benchmark filename, so rule
slugs must never themselves contain a literal `.`.

### Version string

The pipeline version shown in logs is parsed out of the `CHANGELOG` file at Snakemake
startup (`get_version_from_changelog()` in `workflow/Snakefile`), matching a
`Version X.Y.*:` header followed by a leading patch number — bump `CHANGELOG`, not a
separate version constant, when releasing.
