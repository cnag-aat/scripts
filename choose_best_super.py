#!/usr/bin/env python3

# choose_best_super.py
#
# Build the best possible hap1 reference from two haplotype-resolved assemblies
# that are internally phased but NOT interchromosomally phased (typical output of
# hifiasm/Flye + Omni-C curation in PretextViewAI, processed with pretext-to-asm).
#
# For each chromosome (SUPER) the better copy -- judged by a tunable weighted
# composite of BUSCO completeness, contiguity (gaps) and length -- is routed to
# hap1; the worse copy goes to hap2. All user-named sex chromosomes are forced
# into hap1 regardless of score. The hap1<->hap2 correspondence and orientation
# are established from a minimap2 whole-genome alignment (hap2 -> hap1), and hap2
# is renamed/reoriented to match hap1.
#
# Replaces choose_best_super.pl. Reuses the PAF longest-alignment / revcomp logic
# from pafbased_rename.py (Jessica Gomez-Garrido, CNAG).
#
# Author: CNAG Assembly Team
# Contact: tyler.alioto@cnag.eu

import argparse
import os
import re
import shutil
import subprocess
import sys

# --------------------------------------------------------------------------- #
# Small helpers (kept inline so the script is self-contained, mirroring
# pafbased_rename.py).
# --------------------------------------------------------------------------- #

_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'N': 'N', 'n': 'n',
}

_GAP_RE = re.compile(r'[Nn]+')
# SUPER name and optional unloc suffix, e.g. SUPER_3 or SUPER_3_unloc_2
_UNLOC_RE = re.compile(r'^(?P<parent>.+?)_unloc.*$', re.IGNORECASE)
# Trailing integer of a SUPER name, used for canonical ordering / numbering.
_NUM_RE = re.compile(r'(\d+)\s*$')


def revcomp(sequence):
    return ''.join(_COMPLEMENT.get(b, b) for b in reversed(sequence))


def format_sequence(seq, line_length=60):
    return '\n'.join(seq[i:i + line_length] for i in range(0, len(seq), line_length))


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# --------------------------------------------------------------------------- #
# FASTA parsing + per-scaffold stats
# --------------------------------------------------------------------------- #

def read_fasta(path):
    """Yield (name, sequence) tuples. name is the first whitespace-delimited
    token of the header (matches how minimap2 names sequences in a PAF)."""
    if not os.path.isfile(path):
        sys.exit(f"Error: FASTA not found: {path}")
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def gap_stats(seq):
    """Return (gap_count, gap_bases, largest_ungapped_contig) for a sequence."""
    gap_count = 0
    gap_bases = 0
    for m in _GAP_RE.finditer(seq):
        gap_count += 1
        gap_bases += m.end() - m.start()
    # largest stretch of non-N as a proxy for within-SUPER contig contiguity
    largest = 0
    for piece in _GAP_RE.split(seq):
        if len(piece) > largest:
            largest = len(piece)
    return gap_count, gap_bases, largest


def parent_super(name, match):
    """Return the canonical SUPER this scaffold belongs to, or None if it is an
    unplaced scaffold. SUPER_n -> SUPER_n; SUPER_n_unloc_x -> SUPER_n."""
    if match not in name:
        return None
    m = _UNLOC_RE.match(name)
    if m:
        return m.group('parent')
    return name


def super_num(name):
    m = _NUM_RE.search(name)
    return int(m.group(1)) if m else 10 ** 9


def parse_fasta_stats(path, match):
    """Single pass over a haplotype FASTA.

    Returns a dict:
      records[name] = {seq, length, gap_count, gap_bases, largest_contig,
                       parent, is_unloc}
    plus:
      supers[parent] = aggregated stats over the SUPER + its unloc children
      unplaced = [names of scaffolds with no `match` prefix]
    """
    records = {}
    supers = {}
    unplaced = []
    for name, seq in read_fasta(path):
        gc, gb, largest = gap_stats(seq)
        parent = parent_super(name, match)
        is_unloc = parent is not None and parent != name
        records[name] = {
            'seq': seq,
            'length': len(seq),
            'gap_count': gc,
            'gap_bases': gb,
            'largest_contig': largest,
            'parent': parent,
            'is_unloc': is_unloc,
        }
        if parent is None:
            unplaced.append(name)
            continue
        agg = supers.setdefault(parent, {
            'members': [], 'length': 0, 'gap_count': 0,
            'gap_bases': 0, 'largest_contig': 0, 'busco': 0,
        })
        agg['members'].append(name)
        agg['length'] += len(seq)
        agg['gap_count'] += gc
        agg['gap_bases'] += gb
        agg['largest_contig'] = max(agg['largest_contig'], largest)
    return records, supers, unplaced


# --------------------------------------------------------------------------- #
# minimap2 + correspondence
# --------------------------------------------------------------------------- #

def run_minimap2(hap1, hap2, preset, threads, out_paf, dry_run):
    """Align hap2 (query) -> hap1 (target). Returns the command list."""
    cmd = ['minimap2', '-c', '-x', preset, '-t', str(threads), hap1, hap2]
    if dry_run:
        eprint('[dry-run] ' + ' '.join(cmd) + f' > {out_paf}')
        return cmd
    if shutil.which('minimap2') is None:
        sys.exit("Error: minimap2 not found on PATH.")
    eprint('[minimap2] ' + ' '.join(cmd) + f' > {out_paf}')
    with open(out_paf, 'w') as out:
        rc = subprocess.run(cmd, stdout=out).returncode
    if rc != 0:
        sys.exit(f"Error: minimap2 exited with status {rc}")
    return cmd


def correspond_supers(paf_file, match, min_mapq=60):
    """Adapt pafbased_rename.py longest-alignment logic.

    For each hap2 SUPER (query), find the hap1 SUPER (target) holding its single
    longest high-MAPQ aligned block, and the strand of that block.

    Returns:
      h2_to_h1[h2_super] = {'h1': h1_super, 'strand': '+/-', 'len': aln_len}
      h1_claims[h1_super] = [h2_supers...]  (to detect collisions)
    """
    if not os.path.isfile(paf_file):
        sys.exit(f"Error: PAF not found: {paf_file}")
    best = {}  # h2_super -> (h1_super, aln_len, strand)
    with open(paf_file) as fh:
        for ln, line in enumerate(fh, 1):
            hit = line.rstrip('\n').split('\t')
            if len(hit) < 12:
                continue
            q, strand, t = hit[0], hit[4], hit[5]
            try:
                tstart, tend, mapq = int(hit[7]), int(hit[8]), int(hit[11])
            except ValueError:
                continue
            if mapq < min_mapq:
                continue
            qp = parent_super(q, match)
            tp = parent_super(t, match)
            if qp is None or tp is None:
                continue
            # aggregate at the SUPER level using the longest single block
            aln_len = tend - tstart
            if qp not in best or aln_len > best[qp][1]:
                best[qp] = (tp, aln_len, strand)
    h2_to_h1 = {}
    h1_claims = {}
    for h2, (h1, aln_len, strand) in best.items():
        h2_to_h1[h2] = {'h1': h1, 'strand': strand, 'len': aln_len}
        h1_claims.setdefault(h1, []).append(h2)
    return h2_to_h1, h1_claims


# --------------------------------------------------------------------------- #
# BUSCO
# --------------------------------------------------------------------------- #

def run_busco(fasta, lineage, threads, out_path, run_name, dry_run):
    cmd = ['busco', '--miniprot', '-m', 'genome', '-i', fasta,
           '-l', lineage, '-c', str(threads),
           '-o', run_name, '--out_path', out_path, '-f']
    run_dir = os.path.join(out_path, run_name)
    if dry_run:
        eprint('[dry-run] ' + ' '.join(cmd))
        return run_dir
    if shutil.which('busco') is None:
        sys.exit("Error: busco not found on PATH (use --no-busco to skip).")
    eprint('[busco] ' + ' '.join(cmd))
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        sys.exit(f"Error: busco exited with status {rc}")
    return run_dir


def find_full_table(run_dir):
    """Locate full_table.tsv inside a BUSCO run dir (run_<lineage>/full_table.tsv)."""
    direct = os.path.join(run_dir, 'full_table.tsv')
    if os.path.isfile(direct):
        return direct
    for root, _dirs, files in os.walk(run_dir):
        if 'full_table.tsv' in files:
            return os.path.join(root, 'full_table.tsv')
    return None


def busco_per_super(run_dir, match):
    """Parse full_table.tsv -> count Complete (incl. Duplicated) BUSCOs per SUPER.

    full_table.tsv columns: Busco id, Status, Sequence, gene_start, gene_end, ...
    The Sequence field is the scaffold name (may carry a :start-end suffix in
    some BUSCO/miniprot versions, which we strip)."""
    table = find_full_table(run_dir)
    if table is None:
        eprint(f"Warning: no full_table.tsv under {run_dir}; treating BUSCO as 0.")
        return {}
    counts = {}
    with open(table) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            status, seq = f[1], f[2]
            if status not in ('Complete', 'Duplicated'):
                continue
            seq = seq.split(':')[0]
            parent = parent_super(seq, match)
            if parent is None:
                continue
            counts[parent] = counts.get(parent, 0) + 1
    return counts


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #

def _norm(a, b):
    """Min-max normalize two values to [0,1] within the pair. Equal -> (1,1)."""
    if a == b:
        return 1.0, 1.0
    lo, hi = min(a, b), max(a, b)
    return (a - lo) / (hi - lo), (b - lo) / (hi - lo)


def score_pair(s1, s2, weights):
    """Score hap1 SUPER stats `s1` against hap2 SUPER stats `s2`.

    Each metric is normalized within the pair, then combined. Higher = better.
    Contiguity uses negative gap_count (fewer gaps better), tie-broken into the
    largest-ungapped-contig term via the length component already capturing size.
    Returns (score1, score2, components) where components is a dict for the report.
    """
    wb, wc, wl = weights
    # BUSCO (higher better)
    b1n, b2n = _norm(s1['busco'], s2['busco'])
    # contiguity (fewer gaps better -> negate); fall back to largest contig if
    # gap counts tie.
    if s1['gap_count'] == s2['gap_count']:
        c1n, c2n = _norm(s1['largest_contig'], s2['largest_contig'])
    else:
        c1n, c2n = _norm(-s1['gap_count'], -s2['gap_count'])
    # length (longer better)
    l1n, l2n = _norm(s1['length'], s2['length'])
    score1 = wb * b1n + wc * c1n + wl * l1n
    score2 = wb * b2n + wc * c2n + wl * l2n
    comp = {
        'busco_n': (b1n, b2n),
        'contig_n': (c1n, c2n),
        'len_n': (l1n, l2n),
    }
    return score1, score2, comp


# --------------------------------------------------------------------------- #
# Output writing
# --------------------------------------------------------------------------- #

def write_record(out_fh, header, seq, rc=False):
    out_fh.write('>' + header + '\n')
    s = revcomp(seq) if rc else seq
    out_fh.write(format_sequence(s) + '\n')


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def common_prefix(a, b):
    base_a = os.path.basename(a)
    base_b = os.path.basename(b)
    i = 0
    while i < min(len(base_a), len(base_b)) and base_a[i] == base_b[i]:
        i += 1
    pref = base_a[:i].rstrip('_.-')
    return pref or 'assembly'


def main():
    p = argparse.ArgumentParser(
        description="Assign the best SUPER of each chromosome to a hap1 reference "
                    "and the worse to hap2, confirming naming/orientation via "
                    "minimap2 and scoring with BUSCO + contiguity + length.")
    p.add_argument('-1', '--hap1', required=True, help="Haplotype 1 scaffolded FASTA")
    p.add_argument('-2', '--hap2', required=True, help="Haplotype 2 scaffolded FASTA")
    p.add_argument('-l', '--lineage', help="BUSCO lineage (required unless --no-busco)")
    p.add_argument('-x', '--sex', nargs='+', default=[],
                   help="SUPER names always routed to hap1 (e.g. SUPER_X SUPER_Y)")
    p.add_argument('-p', '--prefix', help="Output basename (default: common prefix of inputs)")
    p.add_argument('-o', '--outdir', default='.', help="Output directory (default: .)")
    p.add_argument('-t', '--threads', type=int, default=8, help="Threads (default: 8)")
    p.add_argument('--minimap2-preset', default='asm5', help="minimap2 -x preset (default: asm5)")
    p.add_argument('--paf', help="Reuse an existing PAF instead of running minimap2")
    p.add_argument('--busco-hap1', help="Reuse an existing BUSCO run dir for hap1")
    p.add_argument('--busco-hap2', help="Reuse an existing BUSCO run dir for hap2")
    p.add_argument('--no-busco', action='store_true', help="Skip BUSCO; score on contiguity+length only")
    p.add_argument('--match', default='SUPER', help="Prefix marking a chromosome scaffold (default: SUPER)")
    p.add_argument('--w-busco', type=float, default=3.0, dest='w_busco')
    p.add_argument('--w-contig', type=float, default=2.0, dest='w_contig')
    p.add_argument('--w-len', type=float, default=1.0, dest='w_len')
    p.add_argument('--min-mapq', type=int, default=60, help="Min MAPQ for correspondence (default: 60)")
    p.add_argument('--dry-run', action='store_true',
                   help="Print minimap2/BUSCO commands and the assignment table; write no FASTAs")
    args = p.parse_args()

    if not args.no_busco and not args.lineage and not (args.busco_hap1 and args.busco_hap2):
        sys.exit("Error: --lineage is required unless --no-busco or both --busco-hap1/--busco-hap2 are given.")

    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix or common_prefix(args.hap1, args.hap2)
    sex_set = set(args.sex)
    weights = (args.w_busco, args.w_contig, args.w_len)

    # 1. per-haplotype stats
    eprint(f"[parse] {args.hap1}")
    rec1, sup1, unp1 = parse_fasta_stats(args.hap1, args.match)
    eprint(f"[parse] {args.hap2}")
    rec2, sup2, unp2 = parse_fasta_stats(args.hap2, args.match)

    # 2. minimap2 + correspondence
    paf = args.paf or os.path.join(args.outdir, f"{prefix}.hap2_to_hap1.paf")
    if args.paf and os.path.isfile(args.paf):
        eprint(f"[minimap2] reusing existing PAF: {args.paf}")
    else:
        run_minimap2(args.hap1, args.hap2, args.minimap2_preset,
                     args.threads, paf, args.dry_run)
    if args.dry_run and not os.path.isfile(paf):
        h2_to_h1, h1_claims = {}, {}
    else:
        h2_to_h1, h1_claims = correspond_supers(paf, args.match, args.min_mapq)

    # 3. BUSCO per haplotype -> per-SUPER complete counts
    if not args.no_busco:
        run1 = args.busco_hap1
        run2 = args.busco_hap2
        if not run1:
            run1 = run_busco(args.hap1, args.lineage, args.threads,
                             args.outdir, f"{prefix}_hap1_busco", args.dry_run)
        else:
            eprint(f"[busco] reusing hap1 run dir: {run1}")
        if not run2:
            run2 = run_busco(args.hap2, args.lineage, args.threads,
                             args.outdir, f"{prefix}_hap2_busco", args.dry_run)
        else:
            eprint(f"[busco] reusing hap2 run dir: {run2}")
        if not (args.dry_run and not args.busco_hap1):
            for sup, cnt in busco_per_super(run1, args.match).items():
                if sup in sup1:
                    sup1[sup]['busco'] = cnt
        if not (args.dry_run and not args.busco_hap2):
            for sup, cnt in busco_per_super(run2, args.match).items():
                if sup in sup2:
                    sup2[sup]['busco'] = cnt
    else:
        eprint("Warning: --no-busco; scoring on contiguity + length only.")

    # 4 + 5. build per-chromosome rows and decide.
    # Canonical chromosome key = the hap1 SUPER. For hap2 SUPERs we map them onto
    # their hap1 partner via the correspondence; hap2 SUPERs with no partner stay
    # under their own name.
    rows = []
    used_h2 = set()
    # iterate hap1 SUPERs in numeric order
    for h1 in sorted(sup1, key=super_num):
        partners = [h2 for h2 in h2_to_h1
                    if h2_to_h1[h2]['h1'] == h1 and h2 in sup2]
        # choose the single best (longest aln) hap2 partner if several map here
        h2 = None
        strand = '+'
        if partners:
            partners.sort(key=lambda x: h2_to_h1[x]['len'], reverse=True)
            h2 = partners[0]
            strand = h2_to_h1[h2]['strand']
            used_h2.add(h2)
        s1 = sup1[h1]
        s2 = sup2[h2] if h2 else None
        is_sex = h1 in sex_set or (h2 in sex_set if h2 else False)

        if s2 is None:
            # singleton: only present in hap1
            rows.append({
                'chrom': h1, 'h1_src': h1, 'h2_src': '-', 'strand': '+',
                's1': s1, 's2': None, 'score1': None, 'score2': None,
                'comp': None, 'is_sex': is_sex, 'flag': 'hap1_singleton',
                'winner_hap': 1, 'collision': len(h1_claims.get(h1, [])) > 1,
            })
            continue

        score1, score2, comp = score_pair(s1, s2, weights)
        if is_sex:
            winner = 1  # forced: hap1 copy stays hap1, hap2 partner -> hap2
            flag = 'sex_forced_hap1'
        else:
            winner = 1 if score1 >= score2 else 2
            flag = ''
        rows.append({
            'chrom': h1, 'h1_src': h1, 'h2_src': h2, 'strand': strand,
            's1': s1, 's2': s2, 'score1': score1, 'score2': score2,
            'comp': comp, 'is_sex': is_sex, 'flag': flag,
            'winner_hap': winner,
            'collision': len(h1_claims.get(h1, [])) > 1,
        })

    # hap2 SUPERs that never matched a hap1 SUPER -> singletons on hap2 side
    for h2 in sorted(sup2, key=super_num):
        if h2 in used_h2:
            continue
        is_sex = h2 in sex_set
        rows.append({
            'chrom': h2, 'h1_src': '-', 'h2_src': h2, 'strand': '+',
            's1': None, 's2': sup2[h2], 'score1': None, 'score2': None,
            'comp': None, 'is_sex': is_sex,
            'flag': 'sex_forced_hap1' if is_sex else 'hap2_singleton',
            'winner_hap': 1 if is_sex else 2,
        })

    # ---- assignment report ----
    report_path = os.path.join(args.outdir, f"{prefix}.assignment.tsv")
    corr_path = os.path.join(args.outdir, f"{prefix}.correspondence.tsv")
    header = ['chrom', 'h1_src', 'h2_src', 'strand', 'is_sex',
              'len_h1', 'len_h2', 'gaps_h1', 'gaps_h2', 'busco_h1', 'busco_h2',
              'score_h1', 'score_h2', 'hap1_gets', 'hap2_gets', 'flag']

    def fmt(v):
        if v is None:
            return '-'
        if isinstance(v, float):
            return f"{v:.3f}"
        return str(v)

    report_lines = ['\t'.join(header)]
    for r in rows:
        s1, s2 = r['s1'], r['s2']
        if r['winner_hap'] == 1:
            hap1_gets = r['h1_src'] if r['h1_src'] != '-' else r['h2_src']
            hap2_gets = r['h2_src'] if (s2 and r['h1_src'] != '-') else '-'
        else:
            hap1_gets = r['h2_src']
            hap2_gets = r['h1_src']
        report_lines.append('\t'.join(fmt(x) for x in [
            r['chrom'], r['h1_src'], r['h2_src'], r['strand'], r['is_sex'],
            s1['length'] if s1 else None, s2['length'] if s2 else None,
            s1['gap_count'] if s1 else None, s2['gap_count'] if s2 else None,
            s1['busco'] if s1 else None, s2['busco'] if s2 else None,
            r['score1'], r['score2'], hap1_gets, hap2_gets, r['flag'] or '.',
        ]))

    if args.dry_run:
        eprint("\n[dry-run] assignment table:")
        print('\n'.join(report_lines))
        eprint("\n[dry-run] no FASTA / report files written.")
        return

    with open(report_path, 'w') as fh:
        fh.write('\n'.join(report_lines) + '\n')
    with open(corr_path, 'w') as fh:
        for h2, info in sorted(h2_to_h1.items(), key=lambda kv: super_num(kv[1]['h1'])):
            fh.write(f"{info['h1']}\t{h2}\t{info['strand']}\n")
    eprint(f"[write] {report_path}")
    eprint(f"[write] {corr_path}")

    # ---- write final FASTAs ----
    out1_path = os.path.join(args.outdir, f"{prefix}_hap1.fa")
    out2_path = os.path.join(args.outdir, f"{prefix}_hap2.fa")

    def members_in_order(stats):
        # SUPER first, then its unloc children
        return sorted(stats['members'],
                      key=lambda n: (1 if '_unloc' in n.lower() else 0, n))

    with open(out1_path, 'w') as o1, open(out2_path, 'w') as o2:
        for r in rows:
            canon = r['chrom']  # canonical name = hap1 SUPER (or hap2 singleton name)
            if r['winner_hap'] == 1 and r['h1_src'] != '-':
                # hap1 copy -> hap1 output (no rename/reorient needed)
                for m in members_in_order(sup1[r['h1_src']]):
                    write_record(o1, m, rec1[m]['seq'])
                # hap2 partner -> hap2 output, reoriented + renamed to canon
                if r['h2_src'] != '-' and r['s2'] is not None:
                    rc = (r['strand'] == '-')
                    for m in members_in_order(sup2[r['h2_src']]):
                        new = m.replace(r['h2_src'], canon, 1)
                        write_record(o2, new, rec2[m]['seq'], rc=rc)
            elif r['winner_hap'] == 2:
                # hap2 copy wins -> hap1 output, reoriented + renamed to canon
                rc = (r['strand'] == '-')
                for m in members_in_order(sup2[r['h2_src']]):
                    new = m.replace(r['h2_src'], canon, 1)
                    write_record(o1, new, rec2[m]['seq'], rc=rc)
                # hap1 copy -> hap2 output (keeps its name = canon)
                for m in members_in_order(sup1[r['h1_src']]):
                    write_record(o2, m, rec1[m]['seq'])
            else:
                # hap2-only singleton forced to hap1 (sex) or kept on hap2
                target = o1 if r['winner_hap'] == 1 else o2
                for m in members_in_order(sup2[r['h2_src']]):
                    write_record(target, m, rec2[m]['seq'])

        # unplaced scaffolds: keep with their originating haplotype
        for m in unp1:
            write_record(o1, m, rec1[m]['seq'])
        for m in unp2:
            write_record(o2, m, rec2[m]['seq'])

    eprint(f"[write] {out1_path}")
    eprint(f"[write] {out2_path}")
    eprint("[done]")


if __name__ == '__main__':
    main()
