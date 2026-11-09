#!/usr/bin/env python3
"""
simplify_fasta_headers.py -- short, informative sequence names from ENA/NCBI headers

Reads FASTA (plain or .gz, or '-' for stdin) and writes FASTA with renamed headers:

  >ENA|OZ212012|OZ212012.1 ... chromosome: 1            ->  >CHR_1_OZ212012.1
  >ENA|OZ212015|OZ212015.1 ... chromosome: X            ->  >CHR_X_OZ212015.1
  >ENA|OZ212029|OZ212029.1 ... organelle: mitochondrion ->  >MT_OZ212029.1
  >... organelle: chloroplast / plastid                 ->  >Chloroplast_<acc.ver>
  >ENA|CBAAFG020000022|CBAAFG020000022.1 ... contig: scaffold_68
                                                        ->  >scaffold_68_CBAAFG020000022.1
  >OZ212012.1 ... chromosome 1  (NCBI style)            ->  >CHR_1_OZ212012.1

Headers whose ID is not an accession (pre-publication names such as SUPER_1,
scaffold_12, Chr_3) are left as the first word. Unlocalised/unplaced sequences
are never named as chromosomes.

  simplify_fasta_headers.py asm.fa.gz --map names.tsv > asm.renamed.fa
"""
import argparse
import gzip
import re
import sys

ACC = re.compile(r"^[A-Z]{1,6}_?\d+(\.\d+)?$")


def new_name(header):
    """header: text after '>' (no newline). Returns the new sequence name."""
    parts = header.split(None, 1)
    ident = parts[0]
    desc = parts[1] if len(parts) > 1 else ""
    acc = ident.split("|")[-1]  # ENA|OZ212012|OZ212012.1 -> OZ212012.1
    if not ACC.match(acc):
        return ident  # pre-publication / non-accession name: keep it
    d = desc.lower()

    m = re.search(r"organelle:\s*(\S+)", d) or re.search(r"\b(mitochondri\w*|chloroplast|plastid|apicoplast)\b", d)
    if m:
        org = m.group(1)
        if org.startswith("mito"):
            return f"MT_{acc}"
        if org in ("chloroplast", "plastid"):
            return f"Chloroplast_{acc}"
        return f"{org.capitalize()}_{acc}"

    m = re.search(r"\b(?:contig|scaffold):\s*([^\s,]+)", desc)
    if m:
        return f"{m.group(1)}_{acc}"

    unplaced = re.search(r"unloc|unplaced|unlocali[sz]ed", d)
    m = re.search(r"\bchromosome:?\s+([^\s,]+)", desc, re.I)
    if m and not unplaced:
        return f"CHR_{m.group(1)}_{acc}"
    return acc


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", help="input FASTA (.gz ok), '-' for stdin")
    p.add_argument("--map", help="write new_name<TAB>original_header table here")
    a = p.parse_args()

    if a.fasta == "-":
        fin = sys.stdin
    elif a.fasta.endswith(".gz"):
        fin = gzip.open(a.fasta, "rt")
    else:
        fin = open(a.fasta)
    fmap = open(a.map, "w") if a.map else None
    out = sys.stdout
    seen = set()
    for line in fin:
        if line.startswith(">"):
            hdr = line[1:].rstrip("\n")
            name = new_name(hdr)
            if name in seen:
                sys.exit(f"ERROR: duplicate sequence name '{name}' (from '{hdr}')")
            seen.add(name)
            out.write(f">{name}\n")
            if fmap:
                fmap.write(f"{name}\t{hdr}\n")
        else:
            out.write(line)
    if fmap:
        fmap.close()


if __name__ == "__main__":
    main()
