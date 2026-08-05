#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
from pathlib import Path

def run(cmd, dry_run=False, cwd=None):
    cwd_msg = f" (cwd={cwd})" if cwd else ""
    if dry_run:
        return
    subprocess.run(cmd, check=True, cwd=cwd)

def infer_sample_name(r1_path: Path) -> str:
    name = r1_path.name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    name = re.sub(r'([._-])R?1([._-])', r'\1', name)
    name = re.sub(r'(_R1|\.R1|-R1)$', "", name)
    name = re.sub(r'(?:^|[._-])R1$', "", name)
    name = re.sub(r'(_1)$', "", name)

    return name.rstrip("._-")

def find_pairs(folder: Path, r1_glob: str):
    r1_files = sorted(folder.glob(r1_glob))
    pairs = []
    for r1 in r1_files:
        r1s = r1.name
        candidates = [r1s.replace("_R1", "_R2"),r1s.replace(".R1", ".R2"), r1s.replace("-R1", "-R2"), r1s.replace("R1", "R2"), r1s.replace("_1", "_2"),r1s.replace(".1.", ".2."),]

        r2 = None
        for c in candidates:
            p = folder / c
            if p.exists():
                r2 = p
                break
        if r2 is None:
            print(f"[WARN] No R2 found for R1: {r1.name}")
            continue

        pairs.append((r1, r2))

    return pairs

def main():
    ap = argparse.ArgumentParser(description="Mapping pipeline")
    ap.add_argument("--folder", default=".")
    ap.add_argument("--index", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--count_script", required=True)
    ap.add_argument("--r1_glob", default="*R1*.fastq.gz")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--sort_mem", default="2G")
    ap.add_argument("--outdir", default="bacterial_alignments")
    ap.add_argument("--keep_sorted_bam", action="store_true")
    ap.add_argument("--dry_run", action="store_true")

    args = ap.parse_args()

    folder = Path(args.folder).resolve()
    outdir = (folder / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = Path(args.fasta).resolve()
    count_script = Path(args.count_script).resolve()
    pairs = find_pairs(folder, args.r1_glob)

    if not pairs:
        raise SystemExit("No FASTQ pairs found")

    for r1, r2 in pairs:
        sample = infer_sample_name(r1)
        sam = outdir / f"{sample}_bacterial.sam"
        bam = outdir / f"{sample}_bacterial.bam"
        sbam = outdir / f"{sample}_bacterial.sorted.bam"

        run(["bowtie2", "-x", args.index, "-2", str(r2), "-S", str(sam), "-p", str(args.threads)], args.dry_run)
        run(["samtools", "view", "-bS", str(sam), "-o", str(bam)], args.dry_run)

        if not args.dry_run:
            sam.unlink(missing_ok=True)

        run(["samtools", "sort", "-@", str(args.threads), "-m", args.sort_mem, "-o", str(sbam), str(bam)], args.dry_run)

        if not args.dry_run:
            bam.unlink(missing_ok=True)

        run(["samtools", "index", str(sbam)], args.dry_run)
        run(["python",  str(count_script),  str(fasta), str(sbam)], args.dry_run, cwd=str(outdir))

        if not args.keep_sorted_bam and not args.dry_run:
            bai = Path(str(sbam) + ".bai")
            sbam.unlink(missing_ok=True)
            bai.unlink(missing_ok=True)

        sbam_stem = sbam.stem
        ff=str(fasta)
        print(outdir / f"{sbam_stem}_counts.tsv")
        print(outdir / f"{sbam_stem}_coveredbases.tsv")


if __name__ == "__main__":
    main()