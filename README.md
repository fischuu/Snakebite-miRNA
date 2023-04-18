# Snakebite-miRNA
This is the Luke miRNA pipeline

# Installation
To install, run likt this (using ssh). Please keep in mind, you need to have a registered ssh key at github for this clone to work. If you do not have such a key, pleae proceed via https cloning

```
git clone git@github.com:fischuu/Snakebite-miRNA.git
```

Or via HTTPS

```
git clone https://github.com/fischuu/Snakebite-miRNA.git

```

# Included tools:
 * Bowtie 
 * cutadapt
 * FastQC
 * MultiQC
 * 
# DAG
The visual overview of the various rules of the pipeline
![alt text](https://github.com/fischuu/Pipeline-miRNA/blob/main/workflow.png?raw=true)

# To be implemented
 * Filter the mature.fa and hairpin.fa based on a given species identifier (e.g. bta)
 * Align the reads also to the mature_<identifier>.fa and hairpin_<identifier>.fa 
