#!/usr/bin/env python3
"""
telo_scan.py — Telomeric tandem repeat scanner
================================================
Scans genome assemblies for telomeric motifs with strand separation,
tandem-specificity filtering, softmask awareness, 0-1000 bedgraph
normalisation, and chromosome-end classification.

Outputs (per motif):
  <prefix>.<motif>.5p_telomere.bg      — forward-strand normalised scores
  <prefix>.<motif>.3p_telomere.bg      — reverse-complement normalised scores
  <prefix>.<motif>.5p+3p_telomere.bg   — combined scores
  <prefix>.<motif>.hits.tsv            — every tandem run with coords/copies

Classification outputs (scaffolds >= --min-chr-size):
  <prefix>.classification.tsv          — per-chromosome table + summary stats
  <prefix>.karyogram.strict.png        — graphic, strict classification
  <prefix>.karyogram.loose.png         — graphic, loose classification (orange = loose-only)

Scores are normalised 0-1000 where 1000 = window 100% filled with the motif.
Track names follow PretextGraph convention (5/3 keyboard shortcuts).

Usage:
  python telo_scan.py -i assembly.fasta -m TTAGGG [TTAGGG ...]
      [--window 10000] [--min-tandem 2] [--max-gap 0]
      [--min-chr-size 5000000] [--score-threshold 100] [--terminal-bp 100000]
      [-o prefix] [--threads N] [--no-hits] [--no-karyogram]
"""

import argparse
import gzip
import math
import multiprocessing
import os
import re
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def open_fasta(path: str):
    """Yield (name, sequence) — original case preserved for softmask support."""
    opener = gzip.open if path.endswith(".gz") else open
    name, parts = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def find_tandem_runs(seq: str, motif: str, min_tandem: int, max_gap: int):
    """
    Return list of (start, end, n_copies) for all tandem runs of motif
    that contain >= min_tandem copies, allowing up to max_gap inter-copy bp.
    Matching is case-insensitive; motif must already be uppercase.
    """
    if max_gap == 0:
        pattern = re.compile(
            f"(?:{re.escape(motif)}){{{min_tandem},}}", re.IGNORECASE)
    else:
        spacer = f".{{0,{max_gap}}}"
        pattern = re.compile(
            f"(?:(?:{re.escape(motif)}){spacer}){{{min_tandem - 1},}}{re.escape(motif)}",
            re.IGNORECASE)

    runs, pos = [], 0
    while pos < len(seq):
        m = pattern.search(seq, pos)
        if m is None:
            break
        start, end = m.start(), m.end()
        n_copies = seq[start:end].upper().count(motif)
        if n_copies < min_tandem:
            pos = start + 1
            continue
        runs.append((start, end, n_copies))
        pos = end
    return runs


def build_windows(seq_len: int, window: int):
    for start in range(0, seq_len, window):
        yield start, min(start + window, seq_len)


# ---------------------------------------------------------------------------
# Per-window counting — normalised to 0-1000
# ---------------------------------------------------------------------------

def count_in_windows(runs, seq_len, window, motif_len):
    """
    Sum n_copies per window (apportioning by bp overlap at boundaries) then
    normalise each window to 0-1000, where 1000 = window 100% filled with
    the motif.  The terminal window uses its actual size so it is not
    penalised for being shorter than the nominal window.
    """
    wins = list(build_windows(seq_len, window))
    raw = [0.0] * len(wins)

    for (rstart, rend, ncopies) in runs:
        run_len = rend - rstart
        if run_len == 0:
            continue
        for wi, (ws, we) in enumerate(wins):
            overlap = min(rend, we) - max(rstart, ws)
            if overlap <= 0:
                continue
            raw[wi] += ncopies * overlap / run_len

    scores = []
    for (ws, we), r in zip(wins, raw):
        max_copies = (we - ws) / motif_len
        scores.append(round(r / max_copies * 1000) if max_copies > 0 else 0)

    return wins, scores


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_bedgraph(path, scaffold_windows):
    with open(path, "w") as fh:
        for scaffold, wins, counts in scaffold_windows:
            for (ws, we), cnt in zip(wins, counts):
                fh.write(f"{scaffold}\t{ws}\t{we}\t{cnt}\n")


def write_hits_tsv(path, all_hits):
    with open(path, "w") as fh:
        fh.write("scaffold\tstart\tend\tstrand\tn_copies\tsequence\n")
        for row in all_hits:
            fh.write("\t".join(str(x) for x in row) + "\n")


# ---------------------------------------------------------------------------
# Classification — presence/absence at each chromosome end
# ---------------------------------------------------------------------------

# Karyogram colours
_COL_PRESENT  = "#27AE60"   # green  — telomere confirmed (strict)
_COL_MARGINAL = "#F39C12"   # orange — loose-only (terminal window absent)
_COL_ABSENT   = "#C0392B"   # red    — no signal at this end
_COL_BODY     = "#ECF0F1"   # light-grey chromosome body


def is_sex_chr(name: str) -> bool:
    """Heuristic: True when the scaffold name looks like a sex chromosome."""
    stripped = re.sub(r"^(?:CHR(?:OMOSOME)?_?|LG_?|SUPER_?|SCAFFOLD_?)", "",
                      name.upper())
    return bool(re.match(r"^[XYZW](\d*|_\d+)?$", stripped))


def _terminal_max_score(wins, counts, terminal_bp: int, from_start: bool) -> int:
    """
    Return the highest window score in the terminal region of one end.

    from_start=True  → first terminal_bp of the scaffold (5' / coord start)
    from_start=False → last  terminal_bp of the scaffold (3' / coord end)
    """
    if not wins or not counts:
        return 0
    scaffold_end = wins[-1][1]
    if from_start:
        scores = [c for (ws, _), c in zip(wins, counts) if ws < terminal_bp]
    else:
        scores = [c for (_, we), c in zip(wins, counts)
                  if we > scaffold_end - terminal_bp]
    return max(scores) if scores else 0


def classify_chromosomes(comb_data, scaffold_lengths,
                          min_chr_size, threshold, terminal_bp):
    """
    Classify the 5' (coord start) and 3' (coord end) of every scaffold
    >= min_chr_size.  Uses the combined (fwd+rev) track, taking the
    per-window maximum across all motifs.

    Returns a list of result dicts sorted longest-first.

    Strict: only the outermost single window must exceed threshold.
    Loose:  any window within terminal_bp of the end must exceed threshold.
    """
    # Aggregate: per-window maximum across all motifs
    scaf_wins   = {}
    scaf_counts = {}
    for _motif, entries in comb_data.items():
        for scaffold, wins, counts in entries:
            if scaffold not in scaf_wins:
                scaf_wins[scaffold]   = wins
                scaf_counts[scaffold] = list(counts)
            else:
                scaf_counts[scaffold] = [
                    max(a, b)
                    for a, b in zip(scaf_counts[scaffold], counts)
                ]

    def _pres(score): return score >= threshold

    def _class(s, e):
        if s and e:  return "both"
        if s:        return "5p_only"
        if e:        return "3p_only"
        return "none"

    results = []
    for scaffold, length in scaffold_lengths.items():
        if length < min_chr_size:
            continue
        wins   = scaf_wins.get(scaffold, [])
        counts = scaf_counts.get(scaffold, [])

        ss = counts[0]  if counts else 0   # strict start (first window)
        es = counts[-1] if counts else 0   # strict end   (last window)
        sl = _terminal_max_score(wins, counts, terminal_bp, True)
        el = _terminal_max_score(wins, counts, terminal_bp, False)

        results.append(dict(
            name               = scaffold,
            length             = length,
            length_mb          = length / 1e6,
            # strict
            start_score_strict = ss,
            end_score_strict   = es,
            start_strict       = _pres(ss),
            end_strict         = _pres(es),
            class_strict       = _class(_pres(ss), _pres(es)),
            # loose
            start_score_loose  = sl,
            end_score_loose    = el,
            start_loose        = _pres(sl),
            end_loose          = _pres(el),
            class_loose        = _class(_pres(sl), _pres(el)),
        ))

    results.sort(key=lambda r: -r["length"])
    return results


def write_classification_tsv(path, results, min_chr_size, threshold, terminal_bp):
    """TSV with one row per chromosome and inline summary statistics."""
    n = len(results)

    def pct(k, total):
        return f"{100*k/total:.1f}%" if total else "N/A"

    def summary(mode):
        cats = {"both": 0, "5p_only": 0, "3p_only": 0, "none": 0}
        for r in results:
            cats[r[f"class_{mode}"]] += 1
        return cats

    with open(path, "w") as fh:
        fh.write("# telo_scan telomere classification\n")
        fh.write(f"# min_chr_size={min_chr_size:,} bp  "
                 f"score_threshold={threshold}/1000  "
                 f"terminal_bp={terminal_bp:,}\n")
        fh.write("# 5p = coordinate start of scaffold; "
                 "3p = coordinate end (orient. not guaranteed)\n")
        fh.write("# strict = terminal single window only; "
                 f"loose = any window within {terminal_bp//1000} kb of end\n#\n")

        fh.write("\t".join([
            "scaffold", "length_bp", "length_mb",
            "5p_score_strict", "5p_strict",
            "3p_score_strict", "3p_strict", "class_strict",
            "5p_score_loose",  "5p_loose",
            "3p_score_loose",  "3p_loose",  "class_loose",
        ]) + "\n")

        for r in results:
            fh.write("\t".join([
                r["name"],
                str(r["length"]),
                f"{r['length_mb']:.2f}",
                str(r["start_score_strict"]),
                "Y" if r["start_strict"] else "N",
                str(r["end_score_strict"]),
                "Y" if r["end_strict"]   else "N",
                r["class_strict"],
                str(r["start_score_loose"]),
                "Y" if r["start_loose"]  else "N",
                str(r["end_score_loose"]),
                "Y" if r["end_loose"]    else "N",
                r["class_loose"],
            ]) + "\n")

        for mode, label in [("strict", "terminal window only"),
                             ("loose",  f"terminal {terminal_bp//1000} kb")]:
            c = summary(mode)
            telo = c["both"] + c["5p_only"] + c["3p_only"]
            fh.write(f"#\n# --- SUMMARY ({mode}: {label}) ---\n")
            fh.write(f"# chromosomes analysed (>= {min_chr_size//1_000_000} Mb): {n}\n")
            fh.write(f"# both telomeres : {c['both']:4d}  ({pct(c['both'],   n)})\n")
            fh.write(f"# 5p only        : {c['5p_only']:4d}  ({pct(c['5p_only'], n)})\n")
            fh.write(f"# 3p only        : {c['3p_only']:4d}  ({pct(c['3p_only'], n)})\n")
            fh.write(f"# none           : {c['none']:4d}  ({pct(c['none'],   n)})\n")
            fh.write(f"# >= 1 end       : {telo:4d}  ({pct(telo,         n)})\n")


def draw_karyogram(path, results, mode, title):
    """
    Pill-shaped karyogram.  Chromosomes sorted largest→smallest, autosomes
    first then sex chromosomes.  Each pill cap is coloured by telomere status:

      strict graphic  — green / red
      loose  graphic  — green (strict+loose) / orange (loose-only) / red (absent)
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("[warn] matplotlib not found — skipping karyogram. "
              "Install with: pip install matplotlib", file=sys.stderr)
        return

    autos = sorted([r for r in results if not is_sex_chr(r["name"])],
                   key=lambda r: -r["length"])
    sex   = sorted([r for r in results if is_sex_chr(r["name"])],
                   key=lambda r: -r["length"])
    ordered = autos + sex
    n = len(ordered)
    if n == 0:
        return

    max_mb   = max(r["length_mb"] for r in ordered)
    pill_w   = 0.55
    radius   = pill_w / 2
    x_step   = 1.0
    fig_w    = max(6.0, n * 0.70 + 2.5)
    fig_h    = 9.0

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Use fixed axes bounds rather than a data-unit aspect ratio.  Chromosome
    # length is expressed in Mb, whereas pill width is an ordinal x-axis unit;
    # a circular Wedge in data coordinates would otherwise be flattened almost
    # completely in the rendered image.
    ax.set_xlim(-x_step * 0.9, n * x_step)
    ax.set_ylim(-max_mb * 0.16, max_mb * 1.10)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.18, top=0.93)
    fig.canvas.draw()
    origin_px = ax.transData.transform((0, 0))
    radius_px = ax.transData.transform((radius, 0))[0] - origin_px[0]
    mb_px = ax.transData.transform((0, 1))[1] - origin_px[1]
    cap_h = radius_px / mb_px

    def _cap_colors(r):
        if mode == "strict":
            top = _COL_PRESENT if r["start_strict"] else _COL_ABSENT
            bot = _COL_PRESENT if r["end_strict"]   else _COL_ABSENT
        else:  # loose: 3-way colour
            top = (_COL_PRESENT  if r["start_strict"] else
                   _COL_MARGINAL if r["start_loose"]  else _COL_ABSENT)
            bot = (_COL_PRESENT  if r["end_strict"]   else
                   _COL_MARGINAL if r["end_loose"]    else _COL_ABSENT)
        return top, bot

    for i, r in enumerate(ordered):
        x    = i * x_step
        h    = r["length_mb"]
        ybot = 0.0
        ytop = h
        top_col, bot_col = _cap_colors(r)
        # Keep caps circular in display coordinates, but never let them
        # overlap on very short chromosomes.
        this_cap_h = min(cap_h, h / 2)

        def _arc_points(center_y, start, stop):
            return [
                (x + radius * math.cos(theta),
                 center_y + this_cap_h * math.sin(theta))
                for theta in [start + (stop - start) * j / 24
                              for j in range(25)]
            ]

        top_arc = _arc_points(ytop - this_cap_h, math.pi, 0)
        bot_arc = _arc_points(ybot + this_cap_h, math.pi, 2 * math.pi)

        # Body — no edge so junction with caps is seamless.  The coloured
        # caps occupy the first/last cap height of the chromosome.
        ax.add_patch(mpatches.Rectangle(
            (x - radius, ybot + this_cap_h), pill_w,
            h - 2 * this_cap_h,
            facecolor=_COL_BODY, edgecolor="none", zorder=2))

        # Top cap — 5' / coordinate start
        ax.add_patch(mpatches.Polygon(
            top_arc, closed=True, facecolor=top_col, edgecolor="none",
            zorder=3))

        # Bottom cap — 3' / coordinate end
        ax.add_patch(mpatches.Polygon(
            bot_arc, closed=True, facecolor=bot_col, edgecolor="none",
            zorder=3))

        # Thin outline of the complete pill drawn last (highest zorder)
        ax.plot(*zip(*top_arc), color="#555", linewidth=0.5, zorder=4)
        ax.plot(*zip(*bot_arc), color="#555", linewidth=0.5, zorder=4)
        # Left and right sides
        ax.plot([x - radius, x - radius],
                [ybot + this_cap_h, ytop - this_cap_h],
                color="#555", linewidth=0.5, zorder=4)
        ax.plot([x + radius, x + radius],
                [ybot + this_cap_h, ytop - this_cap_h],
                color="#555", linewidth=0.5, zorder=4)

        # Name label
        ax.text(x, -max_mb * 0.03, r["name"],
                ha="center", va="top", fontsize=12,
                rotation=45, rotation_mode="anchor")

    # Dashed separator between autosomes and sex chromosomes
    # if autos and sex:
    #     sep_x = len(autos) * x_step - x_step * 0.5
    #     ax.axvline(sep_x, color="#95A5A6", linewidth=0.8,
    #                linestyle="--", alpha=0.6, zorder=1)
    #     ax.text(sep_x + 0.12, max_mb * 1.01, "sex",
    #             fontsize=10, color="#7F8C8D", va="bottom")

    # End-type labels on the y-axis side
    ax.text(-0.5, max_mb,  "5′", ha="right", va="center", fontsize=12, color="#444")
    ax.text(-0.5, 0.0,     "3′", ha="right", va="center", fontsize=12, color="#444")

    # Legend
    handles = [
        mpatches.Patch(facecolor=_COL_PRESENT, label="Present (strict)"),
    ]
    if mode == "loose":
        handles.append(
            mpatches.Patch(facecolor=_COL_MARGINAL,
                           label="Present (loose only — terminal window absent)"))
    handles += [
        mpatches.Patch(facecolor=_COL_ABSENT, label="Absent"),
        mpatches.Patch(facecolor=_COL_BODY,   label="Chromosome body"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=12, framealpha=0.92)

    # Axes
    ax.set_ylabel("Length (Mb)", fontsize=12)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v:.0f}"))
    ax.set_xticks([])
    for spine in ("top", "right", "bottom"):
        ax.spines[spine].set_visible(False)

    ax.set_title(title, fontsize=14, pad=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[out] {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Multiprocessing worker
# ---------------------------------------------------------------------------

def _scan_scaffold(args_tuple):
    """Worker: scan one scaffold for all motifs. Returns structured results."""
    scaffold, seq, canonical, min_tandem, max_gap, window, no_hits = args_tuple
    seq_len = len(seq)
    result  = {}
    for fwd_motif, rev_motif in canonical.items():
        motif_len = len(fwd_motif)
        fwd_runs  = find_tandem_runs(seq, fwd_motif, min_tandem, max_gap)
        rev_runs  = find_tandem_runs(seq, rev_motif, min_tandem, max_gap)

        wins_f, cnts_f = count_in_windows(fwd_runs, seq_len, window, motif_len)
        wins_r, cnts_r = count_in_windows(rev_runs, seq_len, window, motif_len)
        cnts_c = [a + b for a, b in zip(cnts_f, cnts_r)]

        hits = []
        if not no_hits:
            for start, end, nc in fwd_runs:
                hits.append((scaffold, start, end, "+", nc, seq[start:end]))
            for start, end, nc in rev_runs:
                hits.append((scaffold, start, end, "-", nc, seq[start:end]))

        result[fwd_motif] = {
            "fwd":  (scaffold, wins_f, cnts_f),
            "rev":  (scaffold, wins_r, cnts_r),
            "comb": (scaffold, wins_f, cnts_c),
            "hits": hits,
        }
    return scaffold, seq_len, result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Telomeric tandem repeat scanner — strand-aware, parallel, "
                    "with end classification and karyogram",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # --- scanning ---
    ap.add_argument("-i", "--input",  required=True,
                    help="Input FASTA (plain or gzipped)")
    ap.add_argument("-m", "--motifs", nargs="+", required=True,
                    help="Telomeric motif(s), 5'→3' forward strand (e.g. TTAGGG)")
    ap.add_argument("--window",     type=int, default=10_000,
                    help="Window size bp (default 10000)")
    ap.add_argument("--min-tandem", type=int, default=2,
                    help="Min consecutive copies to qualify (default 2)")
    ap.add_argument("--max-gap",    type=int, default=0,
                    help="Max extra bases between copies (default 0 = strict tandem)")
    ap.add_argument("-o", "--output", default="telo_out",
                    help="Output prefix (default: telo_out)")
    ap.add_argument("--no-hits",    action="store_true",
                    help="Skip per-run hits TSV")
    ap.add_argument("--threads",    type=int, default=os.cpu_count(),
                    help=f"Worker processes (default: all cores)")

    # --- classification ---
    ap.add_argument("--min-chr-size",    type=int, default=5_000_000,
                    help="Min scaffold length for classification (default: 5,000,000 bp)")
    ap.add_argument("--score-threshold", type=int, default=100,
                    help="Min window score (0-1000) to call telomere present (default: 100)")
    ap.add_argument("--terminal-bp",     type=int, default=100_000,
                    help="Bp from each end for loose classification (default: 100,000)")
    ap.add_argument("--no-karyogram",    action="store_true",
                    help="Skip karyogram PNG output (requires matplotlib)")

    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Canonicalise motifs
    # ------------------------------------------------------------------
    canonical, seen = {}, set()
    for m in [x.upper() for x in args.motifs]:
        rc = revcomp(m)
        if m in seen or rc in seen:
            print(f"[info] Skipping {m} — revcomp already registered",
                  file=sys.stderr)
            continue
        seen.update({m, rc})
        canonical[m] = rc
        print(f"[info] Motif: {m}  revcomp: {rc}", file=sys.stderr)

    fwd_data  = {m: [] for m in canonical}
    rev_data  = {m: [] for m in canonical}
    comb_data = {m: [] for m in canonical}
    hit_data  = {m: [] for m in canonical}

    print(f"[info] input={args.input}  window={args.window}  "
          f"min_tandem={args.min_tandem}  max_gap={args.max_gap}  "
          f"threads={args.threads}", file=sys.stderr)

    def task_iter():
        for scaffold, seq in open_fasta(args.input):
            yield (scaffold, seq, canonical,
                   args.min_tandem, args.max_gap, args.window, args.no_hits)

    scaffold_lengths = {}

    with multiprocessing.Pool(processes=args.threads) as pool:
        for scaffold, seq_len, result in pool.imap(_scan_scaffold, task_iter()):
            print(f"  scaffold {scaffold}  len={seq_len:,}", file=sys.stderr)
            scaffold_lengths[scaffold] = seq_len
            for fwd_motif, data in result.items():
                fwd_data[fwd_motif].append(data["fwd"])
                rev_data[fwd_motif].append(data["rev"])
                comb_data[fwd_motif].append(data["comb"])
                hit_data[fwd_motif].extend(data["hits"])

    # ------------------------------------------------------------------
    # Heuristic: label fwd/rev as 5p or 3p
    # ------------------------------------------------------------------
    TERM_WINS = 5
    MIN_WINS  = 10
    motif_labels = {}
    for fwd_motif in canonical:
        start_score = end_score = 0
        for scaffold, wins, counts in fwd_data[fwd_motif]:
            if len(wins) < MIN_WINS:
                continue
            start_score += sum(counts[:TERM_WINS])
            end_score   += sum(counts[-TERM_WINS:])
        if start_score > end_score:
            fwd_label, rev_label = "5p", "3p"
            reason = f"fwd skews to starts ({start_score} vs {end_score})"
        else:
            fwd_label, rev_label = "3p", "5p"
            reason = (f"fwd skews to ends ({end_score} vs {start_score})"
                      if end_score > start_score
                      else f"tie/no data — using default")
        motif_labels[fwd_motif] = (fwd_label, rev_label)
        print(f"[info] {fwd_motif}: fwd→{fwd_label}_telomere  "
              f"rev→{rev_label}_telomere  ({reason})", file=sys.stderr)

    # ------------------------------------------------------------------
    # Write bedgraphs
    # ------------------------------------------------------------------
    prefix = args.output
    for fwd_motif in canonical:
        tag = fwd_motif
        fwd_label, rev_label = motif_labels[fwd_motif]

        fwd_path  = f"{prefix}.{tag}.{fwd_label}_telomere.bg"
        rev_path  = f"{prefix}.{tag}.{rev_label}_telomere.bg"
        comb_path = f"{prefix}.{tag}.5p+3p_telomere.bg"
        hits_path = f"{prefix}.{tag}.hits.tsv"

        write_bedgraph(fwd_path,  fwd_data[fwd_motif])
        write_bedgraph(rev_path,  rev_data[fwd_motif])
        write_bedgraph(comb_path, comb_data[fwd_motif])
        print(f"[out] {fwd_path}",  file=sys.stderr)
        print(f"[out] {rev_path}",  file=sys.stderr)
        print(f"[out] {comb_path}", file=sys.stderr)

        if not args.no_hits:
            hit_data[fwd_motif].sort(key=lambda x: (x[0], x[1]))
            write_hits_tsv(hits_path, hit_data[fwd_motif])
            print(f"[out] {hits_path}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Classification
    # ------------------------------------------------------------------
    print(f"[info] Classifying ends  min_chr_size={args.min_chr_size:,}  "
          f"threshold={args.score_threshold}  terminal_bp={args.terminal_bp:,}",
          file=sys.stderr)

    cls_results = classify_chromosomes(
        comb_data, scaffold_lengths,
        args.min_chr_size, args.score_threshold, args.terminal_bp,
    )

    if not cls_results:
        print(f"[info] No scaffolds >= {args.min_chr_size:,} bp — "
              "skipping classification outputs.", file=sys.stderr)
    else:
        tsv_path = f"{prefix}.classification.tsv"
        write_classification_tsv(
            tsv_path, cls_results,
            args.min_chr_size, args.score_threshold, args.terminal_bp)
        print(f"[out] {tsv_path}", file=sys.stderr)

        if not args.no_karyogram:
            for mode in ("strict", "loose"):
                if mode == "strict":
                    detail = "terminal window only"
                else:
                    detail = f"terminal {args.terminal_bp // 1000} kb"
                title = (f"Telomere classification — {mode}  "
                         f"(threshold={args.score_threshold}/1000,  {detail})")
                draw_karyogram(
                    f"{prefix}.karyogram.{mode}.png",
                    cls_results, mode, title)

    print("[done]", file=sys.stderr)


if __name__ == "__main__":
    main()
