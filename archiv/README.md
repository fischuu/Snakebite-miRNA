# Snakebite-miRNA
This is the Luke miRNA pipeline

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
