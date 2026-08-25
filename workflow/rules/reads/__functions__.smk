def _reads_sample_row(wildcards):
    """Look up the samples.tsv row matching a sample_id/library_id wildcard pair"""
    matches = samples.query(
        "sample_id == @wildcards.sample_id and library_id == @wildcards.library_id"
    )
    return matches.iloc[0]

def get_raw_read1(wildcards):
    """Path to the raw read1 FASTQ listed in samples.tsv for this sample_id/library_id"""
    return _reads_sample_row(wildcards)["read1"]

def get_concatenation_inputs(wildcards):
    """Every trimmed lane belonging to one sample_id, to be concatenated together"""
    rows = samples.query("sample_id == @wildcards.sample_id")
    return [
        str(TRIMMED / f"{row.sample_id}.{row.library_id}_trimmed.fastq.gz")
        for row in rows.itertuples()
    ]
