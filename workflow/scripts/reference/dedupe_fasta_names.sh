#!/usr/bin/env bash
# Disambiguate duplicate FASTA header names. miRBase occasionally reuses the
# same mature/hairpin name across two paralogous precursor loci within the
# same release, which bowtie-build/STAR happily index but samtools/seqkit
# then refuse to read back ("duplicate reference name" / "duplicate entry
# ... in sam header"). Reads a FASTA on stdin, writes it to stdout: the
# first occurrence of a name is left untouched; the 2nd, 3rd, ... occurrence
# gets "_2", "_3", ... appended to the name itself (not the rest of the
# description), so every reference name is unique before any index is built
# from it.
awk '
/^>/ {
    name = $0
    sub(/^>/, "", name)
    split(name, parts, /[ \t]/)
    id = parts[1]
    count[id]++
    if (count[id] > 1) {
        sub(/^>[^ \t]+/, ">" id "_" count[id])
    }
}
{ print }
'
