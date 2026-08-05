import glob
import sys
import os
import pysam
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

TRIM = 200  
seqs = {}
lengths = {}

for record in SeqIO.parse(sys.argv[1], "fasta"):
    seqs[record.id] = record.description.replace(" ", "_")
    lengths[record.id] = len(record.seq)

results = defaultdict(list)
sample_name = None

for fn in glob.glob(sys.argv[2]):
    print(fn)
    samfile = pysam.AlignmentFile(fn, "rb")

    sample_name = os.path.splitext(os.path.basename(fn))[0]
    print(sample_name)

    for seq in seqs:
        read_count = 0
        covered_positions = set()
        reads = set()
        ref_len = lengths[seq]
        start = TRIM
        end = ref_len - TRIM

        if end <= start:
            continue
        for read in samfile.fetch(seq, start, end):
            ref_positions = read.get_reference_positions()
            read_len = len(ref_positions)

            if read_len == 0:
                continue

            nm = int(read.get_tag("NM")) if read.has_tag("NM") else 0
            read_pid = 1 - (nm / read_len)

            if read_pid >= 0.95 and read_len >= 50 and read.mapping_quality >= 20:
                covered_positions.update(ref_positions)

                if read.query_name not in reads:
                    read_count += 1
                    reads.add(read.query_name)

        results["Sample_Name"].append(sample_name)
        results["reference"].append(seq)
        results["read_count"].append(read_count)
        results["covered_bases"].append(len(covered_positions))

    samfile.close()

results = pd.DataFrame(results)
counts = results.pivot(columns="reference", index="Sample_Name", values="read_count")
print(counts.shape)
counts = counts[counts.columns[counts.sum() > 0]]
print(counts.shape)
counts = counts.reset_index()
counts_file = f"{sample_name}_counts.tsv"
counts.to_csv(counts_file, sep="\t", index=None)
breadth = results.pivot(columns="reference", index="Sample_Name", values="covered_bases")
print(breadth.shape)
breadth = breadth[breadth.columns[breadth.sum() > 0]]
print(breadth.shape)
breadth = breadth.reset_index()
breadth_file = f"{sample_name}_coverage.tsv"
breadth.to_csv(breadth_file, sep="\t", index=None)