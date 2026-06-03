rule reads__link_run:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
    output:
        forward_=READS / "{sample}.{library}.fq.gz",
    log:
        READS / "{sample}.{library}.log"
    benchmark:
        READS / "benchmark/{sample}.{library}.tsv"
    container:
        docker["reads"]
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2>  {log} 1>&2
        """

rule reads__link:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}.fq.gz"
            for sample, library in SAMPLES_LIBRARY
        ],
