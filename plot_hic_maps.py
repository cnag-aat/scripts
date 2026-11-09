#!/usr/bin/env python3
"""
plot_hic_maps.py -- publication-quality Hi-C contact map images from an .mcool

Produces
  * whole-genome maps at several zoom levels (the finest resolution in the
    .mcool that keeps the map under N bins, for each N in --genome-bins)
  * one map per chromosome (finest resolution under --chrom-max-bins bins)

Chromosomes are recognised by name and get a label and boundary lines:
  CHR_1_OZ212012.1 -> "1 (OZ212012.1)"   (as written by simplify_fasta_headers.py)
  SUPER_1, CHR_1, Chr_1, chr_1, chr1     -> "1"   (pre-publication names)
Unlocalised/unplaced sequences (SUPER_12_unloc_1, scaffold_68_..., MT_...) are
drawn but not labelled, or dropped entirely with --chroms-only. If no sequence
is named like a chromosome, a warning is printed and sequences >= --min-size
bp are labelled instead (with their full names). No metadata is needed.

Example
  plot_hic_maps.py asm.mcool -o plots --formats png,pdf
"""
import argparse
import math
import os
import re
import sys

import numpy as np
import cooler
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402

# Juicebox-like white -> red
JUICEBOX = LinearSegmentedColormap.from_list("juicebox", ["#ffffff", "#ff0000"])
# white -> yellow -> red -> black ("fall"-like, cooltools style)
FALL = LinearSegmentedColormap.from_list(
    "fall", ["#ffffff", "#ffeeaa", "#ffcc66", "#ff9933", "#ff3300", "#bb0000", "#330000"])
CMAPS = {"juicebox": JUICEBOX, "fall": FALL}

# CHR_1_OZ212012.1 | CHR_X | SUPER_1 | Chr_3 | chr12 ; not SUPER_12_unloc_1
CHROM_RE = re.compile(r"^(?:CHR|Chr|chr|SUPER)_?([A-Za-z0-9]+)(?:_([A-Z]{1,6}_?\d+\.\d+))?$")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mcool", help="multi-resolution cooler (.mcool)")
    p.add_argument("-o", "--outdir", default="plots", help="output directory [plots]")
    p.add_argument("-p", "--prefix", default=None,
                   help="output file prefix [mcool basename]")
    p.add_argument("--min-size", type=float, default=1e6,
                   help="fallback if no chromosome-named sequences: label seqs >= this many bp [1e6]")
    p.add_argument("--chroms-only", action="store_true",
                   help="drop unlabelled (non-chromosome) sequences from whole-genome maps")
    p.add_argument("--genome-bins", default="1000,2500,4000",
                   help="whole-genome zoom levels as max bins per side [1000,2500,4000]")
    p.add_argument("--genome-resolutions", default=None,
                   help="explicit whole-genome resolutions (bp, comma-sep); overrides --genome-bins")
    p.add_argument("--chrom-max-bins", type=int, default=2500,
                   help="per-chromosome maps: finest resolution with <= this many bins [2500]")
    p.add_argument("--no-chrom-plots", action="store_true", help="skip per-chromosome maps")
    p.add_argument("--raw", action="store_true",
                   help="plot raw counts instead of balanced (ICE) values")
    p.add_argument("--cmap", default="juicebox",
                   help="juicebox, fall, or any matplotlib colormap name [juicebox]")
    p.add_argument("--vmin-pct", type=float, default=5.0,
                   help="lower colour limit as percentile of non-zero pixels [5]")
    p.add_argument("--vmax-pct", type=float, default=99.5,
                   help="upper colour limit as percentile of non-zero pixels [99.5]")
    p.add_argument("--size", type=float, default=10.0, help="figure width/height in inches [10]")
    p.add_argument("--dpi", type=int, default=600, help="raster resolution [600]")
    p.add_argument("--formats", default="png,pdf", help="output formats, comma-sep [png,pdf]")
    p.add_argument("--fontsize", type=float, default=None,
                   help="label font size (default: scaled to the number of chromosomes)")
    p.add_argument("--chrom-regex", default=CHROM_RE.pattern,
                   help="regex identifying chromosome names; group 1 = label, optional "
                        "group 2 = accession shown in parentheses [%(default)s]")
    p.add_argument("--title", default=None, help="title for whole-genome maps")
    return p.parse_args()


def fmt_res(bp):
    if bp >= 1_000_000 and bp % 1_000_000 == 0:
        return f"{bp // 1_000_000}Mb"
    if bp >= 1000 and bp % 1000 == 0:
        return f"{bp // 1000}kb"
    return f"{bp}bp"


def available_resolutions(mcool):
    res = []
    for path in cooler.fileops.list_coolers(mcool):
        m = re.match(r"^/resolutions/(\d+)$", path)
        if m:
            res.append(int(m.group(1)))
    if not res:  # single-resolution .cool
        res = [cooler.Cooler(mcool).binsize]
    return sorted(res)


def open_cooler(mcool, res):
    uri = f"{mcool}::/resolutions/{res}"
    try:
        return cooler.Cooler(uri)
    except Exception:
        return cooler.Cooler(mcool)


def n_bins(lengths, res):
    return sum(math.ceil(L / res) for L in lengths)


def genome_bins(lengths, res):
    """Map width in bins when drawn in genomic coordinates (no per-sequence rounding)."""
    return math.ceil(sum(lengths) / res)


def pick_resolution(resolutions, lengths, max_bins, count=n_bins):
    """Finest resolution whose map has <= max_bins bins; coarsest if none fits."""
    for r in resolutions:
        if count(lengths, r) <= max_bins:
            return r
    return resolutions[-1]


def use_balance(clr, raw):
    if raw:
        return False
    if "weight" in clr.bins().columns:
        return True
    print(f"  [warn] {clr.uri}: no balancing weights, plotting raw counts", file=sys.stderr)
    return False


def colour_limits(mat, vmin_pct, vmax_pct):
    nz = mat[np.isfinite(mat) & (mat > 0)]
    if nz.size == 0:
        return 1e-9, 1.0
    vmin, vmax = np.percentile(nz, [vmin_pct, vmax_pct])
    if vmax <= vmin:
        vmax = nz.max()
        vmin = min(vmin, vmax / 10)
    return float(vmin), float(vmax)


def get_cmap(name):
    cmap = CMAPS.get(name) or matplotlib.colormaps[name]
    cmap = cmap.copy()
    cmap.set_bad("white")
    cmap.set_under("white")
    return cmap


def find_chromosomes(chromsizes, chrom_regex, min_size):
    """Ordered {seq name: label} for the sequences to be labelled."""
    rx = re.compile(chrom_regex)
    labels = {}
    for name in chromsizes.index:
        m = rx.match(name)
        if m:
            acc = m.group(2) if (m.re.groups >= 2) else None
            labels[name] = f"{m.group(1)} ({acc})" if acc else m.group(1)
    if labels:
        print(f"{len(chromsizes)} sequences, {len(labels)} chromosomes (labelled by name)")
        return labels
    print(f"  [warn] no sequence named like a chromosome (/{chrom_regex}/); is this a "
          f"chromosome-level assembly? Labelling sequences >= {min_size:,.0f} bp instead.",
          file=sys.stderr)
    labels = {n: n for n, L in chromsizes.items() if L >= min_size}
    if not labels:
        print("  [warn] no sequence reaches --min-size either: nothing will be labelled",
              file=sys.stderr)
    print(f"{len(chromsizes)} sequences, {len(labels)} >= {min_size:,.0f} bp (labelled)")
    return labels


def draw_map(mat, extent, cmap, vmin, vmax, size, cbar_label):
    fig = plt.figure(figsize=(size, size))
    # [left, bottom, width, height]; room on top/left for labels, right for colourbar
    ax = fig.add_axes([0.10, 0.04, 0.80, 0.80])
    cax = fig.add_axes([0.915, 0.04, 0.015, 0.25])
    im = ax.imshow(mat, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                   interpolation="none", origin="upper", extent=extent, aspect="equal")
    im.set_rasterized(True)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(cbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return fig, ax


def save(fig, base, formats, dpi):
    for fmt in formats:
        out = f"{base}.{fmt}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.1)
        print(f"  wrote {out}")
    plt.close(fig)


def genome_maps(args, resolutions, formats, cmap, chrom_labels):
    clr0 = open_cooler(args.mcool, resolutions[0])
    chromsizes = clr0.chromsizes
    names = list(chromsizes.index)
    keep = [c for c in names if c in chrom_labels or not args.chroms_only]
    keep_len = [int(chromsizes[c]) for c in keep]

    if args.genome_resolutions:
        wanted = [int(float(x)) for x in args.genome_resolutions.split(",")]
        missing = [r for r in wanted if r not in resolutions]
        if missing:
            sys.exit(f"resolutions not in {args.mcool}: {missing}; available: {resolutions}")
    else:
        wanted = [pick_resolution(resolutions, keep_len, int(b), count=genome_bins)
                  for b in args.genome_bins.split(",")]
    wanted = sorted(set(wanted), reverse=True)

    for res in wanted:
        clr = open_cooler(args.mcool, res)
        balance = use_balance(clr, args.raw)
        print(f"whole genome @ {fmt_res(res)} ({genome_bins(keep_len, res)} bins, "
              f"{'balanced' if balance else 'raw'})")

        # bin-index mapping: genome bins of kept sequences -> new contiguous index
        ext = {c: clr.extent(c) for c in names}
        sel = np.concatenate([np.arange(*ext[c]) for c in keep])
        newidx = np.full(clr.info["nbins"], -1, dtype=np.int64)
        newidx[sel] = np.arange(sel.size)

        coo = clr.matrix(balance=balance, sparse=True)[:].tocoo()
        r, c = newidx[coo.row], newidx[coo.col]
        ok = (r >= 0) & (c >= 0) & np.isfinite(coo.data)
        mat = np.zeros((sel.size, sel.size), dtype=np.float32)
        mat[r[ok], c[ok]] = coo.data[ok]
        del coo

        vmin, vmax = colour_limits(mat, args.vmin_pct, args.vmax_pct)

        # Draw in genomic (bp) coordinates, not bins: a sequence shorter than the bin
        # size still has one whole bin, which would inflate small unplaced scaffolds
        # (and stretch chromosomes by up to a bin each). Resample the bin matrix onto
        # a uniform bp pixel grid; each pixel takes the bin under its centre.
        starts, pos = {}, 0
        for ch in keep:
            starts[ch] = pos
            pos += int(chromsizes[ch])
        total = pos
        bins = clr.bins()[:]
        # int64: cooler stores bin ends as int32, and genome offsets exceed 2^31 bp
        # in genomes > 2.1 Gb (int32 arithmetic would silently wrap)
        ends = bins["end"].values.astype(np.int64)
        gen_end = np.concatenate([starts[ch] + ends[ext[ch][0]:ext[ch][1]] for ch in keep])
        if gen_end.size != mat.shape[0] or np.any(np.diff(gen_end) <= 0) or gen_end[-1] != total:
            sys.exit("internal error: genomic bin coordinates are not consistent")
        npix = genome_bins(keep_len, res)
        max_px = int(args.size * 0.8 * args.dpi)  # axes width in output pixels
        npix *= max(1, min(4, max_px // npix))    # oversample so boundaries land finely
        centres = (np.arange(npix) + 0.5) * (total / npix)
        idx = np.minimum(np.searchsorted(gen_end, centres, side="right"), mat.shape[0] - 1)
        mat = mat[np.ix_(idx, idx)]

        fig, ax = draw_map(mat, (0, total, total, 0), cmap, vmin, vmax, args.size,
                           "balanced contacts" if balance else "contacts")

        labelled = [ch for ch in keep if ch in chrom_labels]
        bounds = sorted({starts[ch] for ch in labelled} |
                        {starts[ch] + int(chromsizes[ch]) for ch in labelled})
        lw = 0.3 if len(labelled) > 40 else 0.5
        for b in bounds:
            if 0 < b < total:
                ax.axhline(b, color="0.35", lw=lw)
                ax.axvline(b, color="0.35", lw=lw)

        mids = [starts[ch] + int(chromsizes[ch]) / 2 for ch in labelled]
        labels = [chrom_labels[ch] for ch in labelled]
        fs = args.fontsize or max(3.0, min(9.0, 0.8 * args.size * 72 / max(1, len(labelled)) / 1.6))
        ax.set_xticks(mids)
        ax.set_xticklabels(labels, rotation=90, fontsize=fs)
        ax.set_yticks(mids)
        ax.set_yticklabels(labels, fontsize=fs)
        ax.xaxis.tick_top()
        ax.tick_params(length=0, pad=2)
        ax.set_xlim(0, total)
        ax.set_ylim(total, 0)
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        title = args.title or args.prefix
        fig.text(0.10, 0.005, f"{title}  —  {fmt_res(res)} bins", fontsize=8, va="bottom")

        save(fig, os.path.join(args.outdir, f"{args.prefix}.genome.{fmt_res(res)}"),
             formats, args.dpi)


def chrom_maps(args, resolutions, formats, cmap, chrom_labels):
    clr0 = open_cooler(args.mcool, resolutions[0])
    chroms = list(chrom_labels)
    outdir = os.path.join(args.outdir, "chromosomes")
    os.makedirs(outdir, exist_ok=True)
    for ch in chroms:
        L = int(clr0.chromsizes[ch])
        res = pick_resolution(resolutions, [L], args.chrom_max_bins)
        clr = open_cooler(args.mcool, res)
        balance = use_balance(clr, args.raw)
        mat = clr.matrix(balance=balance).fetch(ch).astype(np.float32)
        mat[~np.isfinite(mat)] = 0
        vmin, vmax = colour_limits(mat, args.vmin_pct, args.vmax_pct)
        mb = L / 1e6
        fig, ax = draw_map(mat, (0, mb, mb, 0), cmap, vmin, vmax, args.size,
                           "balanced contacts" if balance else "contacts")
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.set_xlabel("Mb", fontsize=9)
        ax.set_ylabel("Mb", fontsize=9)
        ax.tick_params(labelsize=8)
        name = chrom_labels[ch]
        fig.text(0.10, 0.005, f"{name}  —  {fmt_res(res)} bins", fontsize=9, va="bottom")
        safe = re.sub(r"[^\w.\-]", "_", ch)
        print(f"{ch} @ {fmt_res(res)}")
        save(fig, os.path.join(outdir, f"{args.prefix}.{safe}.{fmt_res(res)}"),
             formats, args.dpi)


def main():
    args = parse_args()
    if args.prefix is None:
        args.prefix = re.sub(r"\.m?cool$", "", os.path.basename(args.mcool))
    os.makedirs(args.outdir, exist_ok=True)
    formats = [f.strip().lstrip(".") for f in args.formats.split(",") if f.strip()]
    cmap = get_cmap(args.cmap)
    resolutions = available_resolutions(args.mcool)
    print(f"resolutions in {args.mcool}: {', '.join(fmt_res(r) for r in resolutions)}")
    chrom_labels = find_chromosomes(open_cooler(args.mcool, resolutions[0]).chromsizes,
                                    args.chrom_regex, args.min_size)
    genome_maps(args, resolutions, formats, cmap, chrom_labels)
    if not args.no_chrom_plots:
        chrom_maps(args, resolutions, formats, cmap, chrom_labels)


if __name__ == "__main__":
    main()
