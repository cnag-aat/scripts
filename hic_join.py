#!/usr/bin/env python3
"""Apply Hi-C scaffold joins to a chromosome-level assembly, then re-sort and
rename it and regenerate its ENA chromosome / unlocalised lists.

Joins only -- breaks need coordinate-level work and are best done in
PretextView. Each join is written as colon-separated components, each with an
orientation in parentheses ("+" if omitted):

    SUPER_12(-):SUPER_4(+)
    SUPER_3(+):scaffold_47(-):SUPER_9(-)

Joins come from --join (repeatable) and/or --joins FILE (one per line; blank
lines and "#" comments ignored). Components are separated by a run of N
(--gap, default 200).

Classification of input sequences comes from the chromosome and unlocalised
lists, not from names:
  * chromosome-list entries whose label is X/Y/Z/W (optionally numbered: X1,
    X2, ...) are sex chromosomes;
  * entries typed circular-chromosome, located in Mitochondrion/Chloroplast,
    or labelled MT/M/Pltd are organelles;
  * other chromosome-list entries are autosomes;
  * unlocalised-list entries are unlocs of the chromosome they name;
  * everything else in the FASTA is unplaced.

A joined molecule containing a sex chromosome takes that sex label; one
containing an autosome is an autosome; a join of only unplaced scaffolds stays
unplaced (named after its first component). Unlocs of every chromosome in a
join move to the joined chromosome. Any component (unloc or unplaced) used in
a join is consumed.

Output order: autosomes by decreasing length, renamed PREFIX_1..PREFIX_N;
sex chromosomes (X, Y, Z, W; X1 before X2) as PREFIX_<label>; organelles
(names kept); unlocs, grouped by their new chromosome in the order above and
renamed <chrom>_unloc_1..k (k ordered by the old unloc number, components of a
join in join order); unplaced sequences by decreasing length (names kept).

Outputs (given --out-prefix P):
    P.asm.fa.gz                  new assembly (bgzip/pigz if available)
    P.chromosome.list.txt.gz     ENA chromosome list
    P.unlocalized.txt.gz         ENA unlocalised list (if any unlocs)
    P.agp                        AGP 2.1: new sequences in terms of old ones
    P.rename.tsv                 old_name  new_name  orientation  length

Usage:
    hic_join.py -f asm.fa.gz -c chromosome.list.txt.gz -u unlocalized.txt.gz \\
        --join 'SUPER_12(-):SUPER_4(+)' -o ibLobTrog1.3
"""
import argparse
import gzip
import re
import shutil
import subprocess
import sys

SEX_RE = re.compile(r"^(X|Y|Z|W)(\d*)$")
SEX_ORDER = {"X": 0, "Y": 1, "Z": 2, "W": 3}
ORGANELLE_LABELS = {"MT", "M", "MITO", "PLTD", "CP"}
COMP_RE = re.compile(r"^(?P<name>[^()\s]+?)(?:\((?P<strand>[+-])\))?$")
COMP_NUMBERED_UNLOC_RE = re.compile(r"_unloc_(\d+)$")

_COMP = bytes.maketrans(b"ACGTRYKMBVDHNacgtrykmbvdhn",
                        b"TGCAYRMKVBHDNtgcayrmkvbhdn")


def die(msg):
    sys.exit(f"ERROR: {msg}")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """Return ordered dict name -> sequence bytes (case preserved)."""
    if path.endswith(".gz"):
        dec = shutil.which("pigz") or shutil.which("gzip")
        proc = subprocess.Popen([dec, "-dc", path], stdout=subprocess.PIPE,
                                bufsize=1 << 20)
        fh = proc.stdout
    else:
        proc, fh = None, open(path, "rb")
    seqs, name, chunks = {}, None, []
    for line in fh:
        if line.startswith(b">"):
            if name is not None:
                seqs[name] = b"".join(chunks)
            name = line[1:].split()[0].decode()
            if name in seqs:
                die(f"duplicate FASTA name {name}")
            chunks = []
        else:
            chunks.append(line.rstrip(b"\r\n"))
    if name is not None:
        seqs[name] = b"".join(chunks)
    fh.close()
    if proc and proc.wait() != 0:
        die(f"failed to decompress {path}")
    return seqs


def read_chrom_list(path):
    """Return list of [name, label, type, (location...)] rows."""
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3:
                die(f"{path}: expected >=3 tab-separated columns: {line!r}")
            rows.append(f)
    return rows


def read_unloc_list(path):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 2:
                die(f"{path}: expected 2 tab-separated columns: {line!r}")
            rows.append((f[0], f[1]))
    return rows


def parse_join(text):
    comps = []
    for tok in text.strip().split(":"):
        m = COMP_RE.match(tok.strip())
        if not m:
            die(f"cannot parse join component {tok!r} in {text!r}")
        comps.append((m["name"], m["strand"] or "+"))
    if len(comps) < 2:
        die(f"a join needs at least two components: {text!r}")
    return comps


def is_organelle(row):
    label = row[1].upper()
    typ = row[2].lower()
    loc = row[3].lower() if len(row) > 3 else ""
    return (label in ORGANELLE_LABELS or typ.startswith("circular")
            or loc in ("mitochondrion", "chloroplast", "plastid"))


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def unloc_num(name):
    m = COMP_NUMBERED_UNLOC_RE.search(name)
    return int(m.group(1)) if m else 0


def sex_key(label):
    m = SEX_RE.match(label)
    return (SEX_ORDER[m.group(1)], int(m.group(2) or 0))


def open_out_gz(path, threads):
    """Return (writable binary handle, process-or-None)."""
    for tool in ("bgzip", "pigz"):
        exe = shutil.which(tool)
        if exe:
            out = open(path, "wb")
            proc = subprocess.Popen([exe, "-c", "-@" if tool == "bgzip" else "-p",
                                     str(threads)], stdin=subprocess.PIPE,
                                    stdout=out)
            return proc.stdin, (proc, out)
    return gzip.open(path, "wb", compresslevel=6), None


def close_out_gz(fh, handle):
    fh.close()
    if handle:
        proc, out = handle
        if proc.wait() != 0:
            die("compression failed")
        out.close()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See the script header for join syntax and ordering rules.")
    ap.add_argument("-f", "--fasta", required=True, help="assembly FASTA(.gz)")
    ap.add_argument("-c", "--chrom-list", required=True,
                    help="ENA chromosome list (name, label, type[, location])")
    ap.add_argument("-u", "--unloc-list",
                    help="ENA unlocalised list (name, chromosome label)")
    ap.add_argument("-j", "--join", action="append", default=[],
                    help="join, e.g. 'SUPER_12(-):SUPER_4(+)' (repeatable)")
    ap.add_argument("-J", "--joins", help="file with one join per line")
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--gap", type=int, default=200,
                    help="N-gap length between joined components [200]")
    ap.add_argument("--prefix", default="SUPER",
                    help="chromosome name prefix for renaming [SUPER]")
    ap.add_argument("--width", type=int, default=60, help="FASTA line width [60]")
    ap.add_argument("--threads", type=int, default=4,
                    help="compression threads [4]")
    args = ap.parse_args()

    joins_txt = list(args.join)
    if args.joins:
        with open(args.joins) as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    joins_txt.append(line)
    if not joins_txt:
        print("WARNING: no joins given; only sorting/renaming", file=sys.stderr)
    joins = [parse_join(t) for t in joins_txt]

    print(f"Reading {args.fasta} ...", file=sys.stderr)
    seqs = read_fasta(args.fasta)
    chrom_rows = read_chrom_list(args.chrom_list)
    unloc_rows = read_unloc_list(args.unloc_list) if args.unloc_list else []

    # ---- classify input sequences -------------------------------------
    kind = {}          # name -> autosome|sex|organelle|unloc|unplaced
    label_of = {}      # chromosome name -> label
    name_of_label = {}
    organelle_rows = {}
    for row in chrom_rows:
        name, label = row[0], row[1]
        if name not in seqs:
            die(f"chromosome-list entry {name} not in FASTA")
        if label in name_of_label:
            die(f"chromosome label {label} used twice in chromosome list")
        name_of_label[label] = name
        label_of[name] = label
        if is_organelle(row):
            kind[name] = "organelle"
            organelle_rows[name] = row
        elif SEX_RE.match(label):
            kind[name] = "sex"
        else:
            kind[name] = "autosome"
    unloc_parent = {}  # unloc name -> parent chromosome name
    for name, label in unloc_rows:
        if name not in seqs:
            die(f"unlocalised-list entry {name} not in FASTA")
        if label not in name_of_label:
            die(f"unloc {name} refers to chromosome {label} not in chromosome list")
        if name in kind:
            die(f"{name} is in both the chromosome and unlocalised lists")
        kind[name] = "unloc"
        unloc_parent[name] = name_of_label[label]
    for name in seqs:
        kind.setdefault(name, "unplaced")

    # ---- build molecules --------------------------------------------------
    # molecule: dict(components=[(old, strand)], kind, sex_label, unlocs=[old])
    used = {}
    molecules = []
    for comps, txt in zip(joins, joins_txt):
        for name, _ in comps:
            if name not in seqs:
                die(f"join {txt!r}: {name} not in FASTA")
            if name in used:
                die(f"{name} used in two joins ({used[name]!r} and {txt!r})")
            if kind[name] == "organelle":
                die(f"join {txt!r}: cannot join organelle {name}")
            used[name] = txt
        kinds = [kind[n] for n, _ in comps]
        sex_labels = sorted({label_of[n] for n, _ in comps if kind[n] == "sex"})
        if len(sex_labels) > 1:
            die(f"join {txt!r} combines sex chromosomes {sex_labels}")
        if sex_labels:
            mk = "sex"
            if "autosome" in kinds:
                print(f"WARNING: join {txt!r} fuses an autosome onto sex "
                      f"chromosome {sex_labels[0]}; result is labelled "
                      f"{sex_labels[0]}", file=sys.stderr)
        elif "autosome" in kinds:
            mk = "autosome"
        else:
            mk = "unplaced"
            if "unloc" in kinds:
                print(f"WARNING: join {txt!r} has no chromosome component; "
                      "result is unplaced", file=sys.stderr)
        n_chr = kinds.count("autosome") + kinds.count("sex")
        if n_chr > 1:
            print(f"NOTE: join {txt!r} fuses {n_chr} chromosomes",
                  file=sys.stderr)
        molecules.append(dict(components=comps, kind=mk,
                              sex_label=sex_labels[0] if sex_labels else None,
                              join=txt))
    for name in seqs:
        if name in used or kind[name] == "unloc":
            continue
        molecules.append(dict(components=[(name, "+")], kind=kind[name],
                              sex_label=label_of.get(name) if kind[name] == "sex"
                              else None, join=None))
    # unlocs follow their parent chromosome into whatever molecule it is in
    for m in molecules:
        m["unlocs"] = []
        for n, _ in m["components"]:
            if kind[n] in ("autosome", "sex"):
                m["unlocs"].extend(sorted(
                    (u for u, p in unloc_parent.items()
                     if p == n and u not in used), key=unloc_num))
    gap = b"N" * args.gap

    def mol_len(m):
        return (sum(len(seqs[n]) for n, _ in m["components"])
                + args.gap * (len(m["components"]) - 1))

    for m in molecules:
        m["length"] = mol_len(m)

    # ---- order & name -----------------------------------------------------
    auto = sorted((m for m in molecules if m["kind"] == "autosome"),
                  key=lambda m: (-m["length"], m["components"][0][0]))
    sex = sorted((m for m in molecules if m["kind"] == "sex"),
                 key=lambda m: sex_key(m["sex_label"]))
    org = [m for m in molecules if m["kind"] == "organelle"]
    unpl = sorted((m for m in molecules if m["kind"] == "unplaced"),
                  key=lambda m: (-m["length"], m["components"][0][0]))
    for i, m in enumerate(auto, 1):
        m["name"], m["label"] = f"{args.prefix}_{i}", str(i)
    for m in sex:
        m["name"], m["label"] = f"{args.prefix}_{m['sex_label']}", m["sex_label"]
    for m in org:
        m["name"] = m["components"][0][0]
    for m in unpl:
        m["name"] = m["components"][0][0]
    chroms = auto + sex
    unloc_out = []  # (new_name, old_name, parent label)
    for m in chroms:
        for k, u in enumerate(m["unlocs"], 1):
            unloc_out.append((f"{m['name']}_unloc_{k}", u, m["label"]))

    new_names = [m["name"] for m in chroms + org + unpl] + [u[0] for u in unloc_out]
    if len(set(new_names)) != len(new_names):
        die("output name collision (an unplaced sequence already uses a "
            f"{args.prefix}_* name?)")

    # ---- write ------------------------------------------------------------
    P = args.out_prefix
    w = args.width
    fa, fah = open_out_gz(f"{P}.asm.fa.gz", args.threads)
    agp = open(f"{P}.agp", "w")
    agp.write("##agp-version\t2.1\n")
    agp.write(f"# generated by hic_join.py from {args.fasta}\n")
    for t in joins_txt:
        agp.write(f"# join: {t}\n")
    ren = open(f"{P}.rename.tsv", "w")
    ren.write("old_name\tnew_name\torientation\tlength\tjoin\n")

    def emit(new, comps, join_txt):
        parts = []
        pos, part_no = 1, 0
        for i, (n, s) in enumerate(comps):
            if i:
                part_no += 1
                agp.write(f"{new}\t{pos}\t{pos + args.gap - 1}\t{part_no}\tU\t"
                          f"{args.gap}\tscaffold\tyes\tproximity_ligation\n")
                parts.append(gap)
                pos += args.gap
            sq = seqs[n]
            part_no += 1
            agp.write(f"{new}\t{pos}\t{pos + len(sq) - 1}\t{part_no}\tW\t{n}\t1\t"
                      f"{len(sq)}\t{s}\n")
            ren.write(f"{n}\t{new}\t{s}\t{len(sq)}\t{join_txt or ''}\n")
            parts.append(revcomp(sq) if s == "-" else sq)
            pos += len(sq)
        full = b"".join(parts)
        fa.write(f">{new}\n".encode())
        for i in range(0, len(full), w):
            fa.write(full[i:i + w])
            fa.write(b"\n")
        return len(full)

    total_out = 0
    for m in chroms + org:
        total_out += emit(m["name"], m["components"], m["join"])
    for new, old, _ in unloc_out:
        total_out += emit(new, [(old, "+")], None)
    for m in unpl:
        total_out += emit(m["name"], m["components"], m["join"])
    close_out_gz(fa, fah)
    agp.close()
    ren.close()

    total_in = sum(len(s) for s in seqs.values())
    expected = total_in + args.gap * sum(len(j) - 1 for j in joins)
    if total_out != expected:
        die(f"length check failed: wrote {total_out}, expected {expected}")

    with gzip.open(f"{P}.chromosome.list.txt.gz", "wt") as fh:
        for m in chroms:
            fh.write(f"{m['name']}\t{m['label']}\tchromosome\n")
        for m in org:
            fh.write("\t".join(organelle_rows[m["name"]]) + "\n")
    if unloc_out:
        with gzip.open(f"{P}.unlocalized.txt.gz", "wt") as fh:
            for new, _, label in unloc_out:
                fh.write(f"{new}\t{label}\n")

    # ---- report -------------------------------------------------------------
    print(f"Input: {len(seqs)} seqs, {total_in:,} bp; output "
          f"{len(new_names)} seqs, {total_out:,} bp "
          f"(+{expected - total_in} bp of gap in {len(joins)} join(s))",
          file=sys.stderr)
    print(f"{'new':<22}{'length':>14}  from", file=sys.stderr)
    for m in chroms + org:
        src = ":".join(f"{n}({s})" for n, s in m["components"])
        flag = "  <- joined" if m["join"] else ""
        print(f"{m['name']:<22}{m['length']:>14,}  {src}{flag}", file=sys.stderr)
    for new, old, _ in unloc_out:
        if new != old:
            print(f"{new:<22}{len(seqs[old]):>14,}  {old}", file=sys.stderr)
    for m in unpl:
        if m["join"]:
            src = ":".join(f"{n}({s})" for n, s in m["components"])
            print(f"{m['name']:<22}{m['length']:>14,}  {src}  <- joined",
                  file=sys.stderr)
    print(f"{len(unpl)} unplaced sequences (names unchanged, longest first)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
