#!/usr/bin/env python3

# choose_best_super.py
#
# Build the best possible hap1 reference from two haplotype-resolved assemblies
# that are internally phased but NOT interchromosomally phased (typical output of
# hifiasm/Flye + Omni-C curation in PretextViewAI, processed with pretext-to-asm).
#
# For each autosome (SUPER) the better copy -- judged by a tunable weighted
# composite of BUSCO completeness, contiguity (gaps) and length -- is routed to
# hap1; the worse copy goes to hap2. ALL sex chromosomes (SUPER_X*, SUPER_Y*,
# SUPER_Z*, SUPER_W*, plus any named with --sex) are routed to hap1 regardless of
# haplotype of origin. Sex chromosomes are matched across haplotypes by name, NOT
# by alignment, so a Z is never cross-matched/renamed to a W.
#
# Autosome hap1<->hap2 correspondence and orientation come from a minimap2
# whole-genome alignment (hap2 -> hap1). Output sequence names are the bare
# chromosome token (e.g. SUPER_3, SUPER_W, SUPER_3_unloc_1) with the original
# haplotype prefix stripped; the prefix/origin/orientation is recorded in the
# lookup table.
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
# Small helpers
# --------------------------------------------------------------------------- #

_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'N': 'N', 'n': 'n',
}

_GAP_RE = re.compile(r'[Nn]+')
_UNLOC_RE = re.compile(r'^(?P<parent>.+?)_unloc.*$', re.IGNORECASE)
_NUM_RE = re.compile(r'(\d+)\s*$')


def revcomp(sequence):
    return ''.join(_COMPLEMENT.get(b, b) for b in reversed(sequence))


def format_sequence(seq, line_length=60):
    return '\n'.join(seq[i:i + line_length] for i in range(0, len(seq), line_length))


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# --------------------------------------------------------------------------- #
# Chromosome-token helpers
#
# Sequences are named like "rHemHip.H1.SUPER_3" or "rHemHip.H1.SUPER_W_unloc_1".
# The *token* is the part starting at the match prefix ("SUPER_..."); the
# *parent token* drops any _unloc suffix; the *prefix* is everything before the
# token (the per-haplotype label we strip from output and keep in the lookup).
# --------------------------------------------------------------------------- #

def chrom_token(name, match):
    i = name.find(match)
    return name[i:] if i != -1 else None


def parent_token(token):
    if token is None:
        return None
    m = _UNLOC_RE.match(token)
    return m.group('parent') if m else token


def name_prefix(name, match):
    i = name.find(match)
    return name[:i] if i != -1 else ''


def member_suffix(member_token, src_parent):
    """The part of a member token after its parent (e.g. '_unloc_1' or '')."""
    return member_token[len(src_parent):]


def super_num(token):
    m = _NUM_RE.search(token)
    return int(m.group(1)) if m else 10 ** 9


def make_sex_test(match, sex_prefixes, extra_names, auto):
    """Return is_sex(parent_token) -> bool."""
    extra = set(extra_names or [])
    pat = None
    if auto and sex_prefixes:
        letters = ''.join(re.escape(p) for p in sex_prefixes)
        pat = re.compile(rf'^{re.escape(match)}_[{letters}]\d*$', re.IGNORECASE)

    def is_sex(ptoken):
        if ptoken in extra:
            return True
        return bool(pat and pat.match(ptoken))

    return is_sex


# --------------------------------------------------------------------------- #
# FASTA parsing + per-scaffold stats
# --------------------------------------------------------------------------- #

def read_fasta(path):
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
    gap_count = gap_bases = 0
    for m in _GAP_RE.finditer(seq):
        gap_count += 1
        gap_bases += m.end() - m.start()
    largest = max((len(p) for p in _GAP_RE.split(seq)), default=0)
    return gap_count, gap_bases, largest


def parse_fasta_stats(path, match):
    """Single pass over a haplotype FASTA. Keys SUPERs by parent token.

    Returns:
      records[fullname] = {seq, length, token}
      supers[parent_token] = {members:[fullnames], length, gap_count, gap_bases,
                              largest_contig, busco, prefix}
      unplaced = [fullnames without the match prefix]
    """
    records, supers, unplaced = {}, {}, []
    for name, seq in read_fasta(path):
        token = chrom_token(name, match)
        records[name] = {'seq': seq, 'length': len(seq), 'token': token}
        if token is None:
            unplaced.append(name)
            continue
        gc, gb, largest = gap_stats(seq)
        ptoken = parent_token(token)
        agg = supers.setdefault(ptoken, {
            'members': [], 'length': 0, 'gap_count': 0, 'gap_bases': 0,
            'largest_contig': 0, 'busco': 0, 'prefix': name_prefix(name, match),
        })
        agg['members'].append(name)
        agg['length'] += len(seq)
        agg['gap_count'] += gc
        agg['gap_bases'] += gb
        agg['largest_contig'] = max(agg['largest_contig'], largest)
    return records, supers, unplaced


def members_in_order(stats):
    return sorted(stats['members'], key=lambda n: (1 if '_unloc' in n.lower() else 0, n))


# --------------------------------------------------------------------------- #
# minimap2 + correspondence
# --------------------------------------------------------------------------- #

def run_minimap2(hap1, hap2, preset, threads, out_paf, dry_run):
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


def correspond_supers(paf_file, match, is_sex, min_mapq=60):
    """Longest high-MAPQ block per query token -> target token.

    Sex chromosomes are EXCLUDED from autosome correspondence so they are never
    cross-matched (e.g. Z<->W). A separate self_strand map records orientation
    for same-token alignments (used to reorient paired sex chromosomes).

    Returns:
      h2_to_h1[h2_parent] = {'h1': h1_parent, 'strand': '+/-', 'len': aln_len}
      h1_claims[h1_parent] = [h2_parents...]
      self_strand[parent]  = '+/-' for the longest same-token block
    """
    if not os.path.isfile(paf_file):
        sys.exit(f"Error: PAF not found: {paf_file}")
    best, self_best = {}, {}
    with open(paf_file) as fh:
        for line in fh:
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
            qp = parent_token(chrom_token(q, match))
            tp = parent_token(chrom_token(t, match))
            if qp is None or tp is None:
                continue
            aln_len = tend - tstart
            if qp == tp and (qp not in self_best or aln_len > self_best[qp][0]):
                self_best[qp] = (aln_len, strand)
            if is_sex(qp) or is_sex(tp):
                continue
            if qp not in best or aln_len > best[qp][1]:
                best[qp] = (tp, aln_len, strand)
    h2_to_h1, h1_claims = {}, {}
    for h2, (h1, aln_len, strand) in best.items():
        h2_to_h1[h2] = {'h1': h1, 'strand': strand, 'len': aln_len}
        h1_claims.setdefault(h1, []).append(h2)
    self_strand = {p: s for p, (_l, s) in self_best.items()}
    return h2_to_h1, h1_claims, self_strand


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
    direct = os.path.join(run_dir, 'full_table.tsv')
    if os.path.isfile(direct):
        return direct
    for root, _dirs, files in os.walk(run_dir):
        if 'full_table.tsv' in files:
            return os.path.join(root, 'full_table.tsv')
    return None


def busco_per_super(run_dir, match):
    """Parse full_table.tsv -> Complete (incl. Duplicated) BUSCOs per parent token."""
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
            if len(f) < 3 or f[1] not in ('Complete', 'Duplicated'):
                continue
            ptoken = parent_token(chrom_token(f[2].split(':')[0], match))
            if ptoken is None:
                continue
            counts[ptoken] = counts.get(ptoken, 0) + 1
    return counts


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #

def _norm(a, b):
    if a == b:
        return 1.0, 1.0
    lo, hi = min(a, b), max(a, b)
    return (a - lo) / (hi - lo), (b - lo) / (hi - lo)


def score_pair(s1, s2, weights):
    wb, wc, wl = weights
    b1n, b2n = _norm(s1['busco'], s2['busco'])
    if s1['gap_count'] == s2['gap_count']:
        c1n, c2n = _norm(s1['largest_contig'], s2['largest_contig'])
    else:
        c1n, c2n = _norm(-s1['gap_count'], -s2['gap_count'])
    l1n, l2n = _norm(s1['length'], s2['length'])
    return wb * b1n + wc * c1n + wl * l1n, wb * b2n + wc * c2n + wl * l2n


# --------------------------------------------------------------------------- #
# Output writing
# --------------------------------------------------------------------------- #

def write_record(out_fh, header, seq, rc=False):
    out_fh.write('>' + header + '\n')
    out_fh.write(format_sequence(revcomp(seq) if rc else seq) + '\n')


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def common_prefix(a, b):
    ba, bb = os.path.basename(a), os.path.basename(b)
    i = 0
    while i < min(len(ba), len(bb)) and ba[i] == bb[i]:
        i += 1
    pref = ba[:i].rstrip('_.-')
    # strip a trailing haplotype tag like ".H" (rHemHip.H1/.H2 -> rHemHip.H -> rHemHip)
    if pref.endswith('.H'):
        pref = pref[:-2]
    return pref.rstrip('_.-') or 'assembly'


def main():
    p = argparse.ArgumentParser(
        description="Assign the best SUPER of each autosome to a hap1 reference and "
                    "the worse to hap2 (scored by BUSCO + contiguity + length), routing "
                    "ALL sex chromosomes to hap1. Confirms autosome naming/orientation "
                    "via minimap2; sex chromosomes are matched by name, never alignment.")
    p.add_argument('-1', '--hap1', required=True, help="Haplotype 1 scaffolded FASTA")
    p.add_argument('-2', '--hap2', required=True, help="Haplotype 2 scaffolded FASTA")
    p.add_argument('-l', '--lineage', help="BUSCO lineage (required unless --no-busco)")
    p.add_argument('-x', '--sex', nargs='+', default=[],
                   help="Extra SUPER tokens to treat as sex chromosomes (e.g. SUPER_B1). "
                        "Combined with auto-detection of SUPER_X*/Y*/Z*/W*.")
    p.add_argument('--sex-prefixes', nargs='+', default=['X', 'Y', 'Z', 'W'],
                   help="Single-letter chromosome prefixes auto-classified as sex "
                        "(default: X Y Z W). SUPER_<letter><digits?> matches.")
    p.add_argument('--no-auto-sex', action='store_true',
                   help="Disable pattern-based sex detection; use only --sex names.")
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
    weights = (args.w_busco, args.w_contig, args.w_len)
    is_sex = make_sex_test(args.match, args.sex_prefixes, args.sex, not args.no_auto_sex)

    # 1. per-haplotype stats (keyed by parent token)
    eprint(f"[parse] {args.hap1}")
    rec1, sup1, unp1 = parse_fasta_stats(args.hap1, args.match)
    eprint(f"[parse] {args.hap2}")
    rec2, sup2, unp2 = parse_fasta_stats(args.hap2, args.match)

    # 2. minimap2 + autosome correspondence (sex excluded)
    paf = args.paf or os.path.join(args.outdir, f"{prefix}.hap2_to_hap1.paf")
    if args.paf and os.path.isfile(args.paf):
        eprint(f"[minimap2] reusing existing PAF: {args.paf}")
    else:
        run_minimap2(args.hap1, args.hap2, args.minimap2_preset, args.threads, paf, args.dry_run)
    if args.dry_run and not os.path.isfile(paf):
        h2_to_h1, h1_claims, self_strand = {}, {}, {}
    else:
        h2_to_h1, h1_claims, self_strand = correspond_supers(paf, args.match, is_sex, args.min_mapq)

    # 3. BUSCO per haplotype -> per-token complete counts
    if not args.no_busco:
        run1, run2 = args.busco_hap1, args.busco_hap2
        if not run1:
            run1 = run_busco(args.hap1, args.lineage, args.threads, args.outdir,
                             f"{prefix}_hap1_busco", args.dry_run)
        else:
            eprint(f"[busco] reusing hap1 run dir: {run1}")
        if not run2:
            run2 = run_busco(args.hap2, args.lineage, args.threads, args.outdir,
                             f"{prefix}_hap2_busco", args.dry_run)
        else:
            eprint(f"[busco] reusing hap2 run dir: {run2}")
        if not (args.dry_run and not args.busco_hap1):
            for tok, cnt in busco_per_super(run1, args.match).items():
                if tok in sup1:
                    sup1[tok]['busco'] = cnt
        if not (args.dry_run and not args.busco_hap2):
            for tok, cnt in busco_per_super(run2, args.match).items():
                if tok in sup2:
                    sup2[tok]['busco'] = cnt
    else:
        eprint("Warning: --no-busco; scoring on contiguity + length only.")

    # 4 + 5. build per-chromosome rows and decide.
    # winner_hap = which haplotype's copy goes into the hap1 reference.
    rows = []
    used_h2 = set()

    # --- autosomes: hap1 tokens with a minimap2 partner, or hap1-only ---
    for h1 in sorted([t for t in sup1 if not is_sex(t)], key=super_num):
        partners = [h2 for h2 in h2_to_h1
                    if h2_to_h1[h2]['h1'] == h1 and h2 in sup2 and not is_sex(h2)]
        if partners:
            partners.sort(key=lambda x: h2_to_h1[x]['len'], reverse=True)
            h2 = partners[0]
            strand = h2_to_h1[h2]['strand']
            used_h2.add(h2)
            s1, s2 = sup1[h1], sup2[h2]
            sc1, sc2 = score_pair(s1, s2, weights)
            rows.append({
                'chrom': h1, 'kind': 'pair', 'h1': h1, 'h2': h2, 'strand': strand,
                's1': s1, 's2': s2, 'score1': sc1, 'score2': sc2, 'is_sex': False,
                'winner_hap': 1 if sc1 >= sc2 else 2, 'flag': '',
                'collision': len(h1_claims.get(h1, [])) > 1,
            })
        else:
            rows.append({
                'chrom': h1, 'kind': 'h1_single', 'h1': h1, 'h2': None, 'strand': '+',
                's1': sup1[h1], 's2': None, 'score1': None, 'score2': None,
                'is_sex': False, 'winner_hap': 1, 'flag': 'hap1_only', 'collision': False,
            })

    # --- autosomes only in hap2 (no confident hap1 partner) -> hap1 reference ---
    for h2 in sorted([t for t in sup2 if not is_sex(t)], key=super_num):
        if h2 in used_h2:
            continue
        rows.append({
            'chrom': h2, 'kind': 'h2_single', 'h1': None, 'h2': h2, 'strand': '+',
            's1': None, 's2': sup2[h2], 'score1': None, 'score2': None,
            'is_sex': False, 'winner_hap': 1, 'flag': 'hap2_only_to_hap1', 'collision': False,
        })

    # --- sex chromosomes: matched by token name across haps, all -> hap1 ---
    sex_tokens = sorted(
        {t for t in sup1 if is_sex(t)} | {t for t in sup2 if is_sex(t)},
        key=lambda t: (super_num(t), t))
    for tok in sex_tokens:
        in1, in2 = tok in sup1, tok in sup2
        if in1 and in2:
            s1, s2 = sup1[tok], sup2[tok]
            sc1, sc2 = score_pair(s1, s2, weights)
            rows.append({
                'chrom': tok, 'kind': 'pair', 'h1': tok, 'h2': tok,
                'strand': self_strand.get(tok, '+'),
                's1': s1, 's2': s2, 'score1': sc1, 'score2': sc2, 'is_sex': True,
                'winner_hap': 1 if sc1 >= sc2 else 2, 'flag': 'sex_both_haps',
                'collision': False,
            })
        elif in1:
            rows.append({
                'chrom': tok, 'kind': 'h1_single', 'h1': tok, 'h2': None, 'strand': '+',
                's1': sup1[tok], 's2': None, 'score1': None, 'score2': None,
                'is_sex': True, 'winner_hap': 1, 'flag': 'sex_hap1', 'collision': False,
            })
        else:
            rows.append({
                'chrom': tok, 'kind': 'h2_single', 'h1': None, 'h2': tok, 'strand': '+',
                's1': None, 's2': sup2[tok], 'score1': None, 'score2': None,
                'is_sex': True, 'winner_hap': 1, 'flag': 'sex_hap2', 'collision': False,
            })

    # ---- assignment report ----
    header = ['chrom', 'h1_src', 'h2_src', 'strand', 'is_sex',
              'len_h1', 'len_h2', 'gaps_h1', 'gaps_h2', 'busco_h1', 'busco_h2',
              'score_h1', 'score_h2', 'hap1_gets', 'hap2_gets', 'flag']

    def fmt(v):
        if v is None:
            return '-'
        return f"{v:.3f}" if isinstance(v, float) else str(v)

    report_lines = ['\t'.join(header)]
    for r in rows:
        s1, s2 = r['s1'], r['s2']
        if r['kind'] == 'pair':
            hap1_gets, hap2_gets = (r['h1'], r['h2']) if r['winner_hap'] == 1 else (r['h2'], r['h1'])
        elif r['kind'] == 'h1_single':
            hap1_gets, hap2_gets = r['h1'], '-'
        else:
            hap1_gets, hap2_gets = r['h2'], '-'
        report_lines.append('\t'.join(fmt(x) for x in [
            r['chrom'], r['h1'] or '-', r['h2'] or '-', r['strand'], r['is_sex'],
            s1['length'] if s1 else None, s2['length'] if s2 else None,
            s1['gap_count'] if s1 else None, s2['gap_count'] if s2 else None,
            s1['busco'] if s1 else None, s2['busco'] if s2 else None,
            r['score1'], r['score2'], hap1_gets, hap2_gets, r['flag'] or '.',
        ]))

    if args.dry_run:
        eprint("\n[dry-run] assignment table:")
        print('\n'.join(report_lines))
        eprint("\n[dry-run] no FASTA / lookup files written.")
        return

    out1_path = os.path.join(args.outdir, f"{prefix}.reassigned.hap1.fa")
    out2_path = os.path.join(args.outdir, f"{prefix}.reassigned.hap2.fa")
    report_path = os.path.join(args.outdir, f"{prefix}.reassigned.assignment.tsv")
    lookup_path = os.path.join(args.outdir, f"{prefix}.reassigned.lookup.tsv")

    with open(report_path, 'w') as fh:
        fh.write('\n'.join(report_lines) + '\n')

    # ---- write final FASTAs + lookup ----
    lookup = [['new_name', 'orig_name', 'source_hap', 'orientation', 'dest']]
    seen_hap1, seen_hap2 = {}, {}

    def emit(out_fh, dest, canon, src_stats, rec, src_hap, rc):
        """Write all members of src_stats to out_fh, renamed to `canon` parent,
        with rc applied. Records lookup rows and tracks output-name collisions."""
        src_parent = src_stats['members'][0]
        src_parent = parent_token(rec[src_parent]['token'])
        for m in members_in_order(src_stats):
            suffix = member_suffix(rec[m]['token'], src_parent)
            new = canon + suffix
            seen = seen_hap1 if dest == 'hap1' else seen_hap2
            seen[new] = seen.get(new, 0) + 1
            write_record(out_fh, new, rec[m]['seq'], rc=rc)
            lookup.append([new, m, src_hap, '-' if rc else '+', dest])

    with open(out1_path, 'w') as o1, open(out2_path, 'w') as o2:
        for r in rows:
            canon = r['chrom']
            rc = (r['strand'] == '-')
            if r['kind'] == 'pair':
                if r['winner_hap'] == 1:
                    emit(o1, 'hap1', canon, sup1[r['h1']], rec1, 'hap1', False)
                    emit(o2, 'hap2', canon, sup2[r['h2']], rec2, 'hap2', rc)
                else:
                    emit(o1, 'hap1', canon, sup2[r['h2']], rec2, 'hap2', rc)
                    emit(o2, 'hap2', canon, sup1[r['h1']], rec1, 'hap1', False)
            elif r['kind'] == 'h1_single':
                emit(o1, 'hap1', canon, sup1[r['h1']], rec1, 'hap1', False)
            else:  # h2_single -> hap1
                emit(o1, 'hap1', canon, sup2[r['h2']], rec2, 'hap2', False)

        for m in unp1:
            write_record(o1, m, rec1[m]['seq'])
            lookup.append([m, m, 'hap1', '+', 'hap1'])
        for m in unp2:
            write_record(o2, m, rec2[m]['seq'])
            lookup.append([m, m, 'hap2', '+', 'hap2'])

    with open(lookup_path, 'w') as fh:
        for row in lookup:
            fh.write('\t'.join(row) + '\n')

    # collision warnings (duplicate output names within a haplotype)
    dups1 = sorted(n for n, c in seen_hap1.items() if c > 1)
    dups2 = sorted(n for n, c in seen_hap2.items() if c > 1)
    if dups1:
        eprint(f"WARNING: duplicate names in hap1 output (ambiguous correspondence): {dups1}")
    if dups2:
        eprint(f"WARNING: duplicate names in hap2 output: {dups2}")

    eprint(f"[write] {report_path}")
    eprint(f"[write] {lookup_path}")
    eprint(f"[write] {out1_path}")
    eprint(f"[write] {out2_path}")
    eprint("[done]")


if __name__ == '__main__':
    main()
