#!/usr/bin/env python3
import sys, os

bin_table, assembly_fasta, outdir = sys.argv[1:4]
os.makedirs(outdir, exist_ok=True)

contig2bin = {}
with open(bin_table) as f:
    for line in f:
        cols = line.rstrip("\n").split("\t")
        bin_id, contig = cols[0], cols[1]
        bin_id = os.path.basename(bin_id)   # strip any path (e.g. "magout/") baked into the bin name
        contig2bin[contig] = bin_id

handles = {}
def get_handle(bin_id):
    if bin_id not in handles:
        handles[bin_id] = open(os.path.join(outdir, f"{bin_id}.fa"), "w")
    return handles[bin_id]

current_bin = None
with open(assembly_fasta) as fasta:
    for line in fasta:
        if line.startswith(">"):
            contig_id = line[1:].split()[0].strip()
            current_bin = contig2bin.get(contig_id)
            if current_bin:
                get_handle(current_bin).write(line)
        elif current_bin:
            get_handle(current_bin).write(line)

for h in handles.values():
    h.close()
