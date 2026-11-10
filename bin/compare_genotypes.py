#!/usr/bin/env python3
"""
compare_genotypes.py
──────────────────────────────────────────────────────────────────────────────
Cross-sample genotype concordance from cellsnp-lite pseudo-bulk VCFs.

Flags two types of problems:
  ● RED border   — same_patient=True  + concordance LOW  (<0.95)
                   → samples labeled as same donor but genotypes differ
                   → likely swap / mislabeling within a patient
  ● ORANGE border — same_patient=False + concordance HIGH (≥0.99)
                   → samples labeled as different donors but genotypes match
                   → likely swap across patients or contamination

Outputs:
  concordance_all_pairs.tsv         – every pairwise comparison
  concordance_same_patient.tsv      – only within-patient pairs
  concordance_flags.tsv             – flagged pairs only
  genotype_concordance_heatmap.pdf  – publication-quality PDF
"""
import argparse, gzip, re, sys
from pathlib import Path
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import seaborn as sns


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf_dir",       default=None,
                   help="Directory with cellSNP.*.vcf.gz files (required unless --skip_genotyping)")
    p.add_argument("--samplesheet",   required=True)
    p.add_argument("--skip_genotyping", action='store_true',
                   help="Skip VCF parsing and concordance; reload concordance_all_pairs.tsv "
                        "from --outdir and regenerate plots only")
    p.add_argument("--min_dp",        type=int,   default=20)
    p.add_argument("--max_oth_frac",  type=float, default=0.10)
    p.add_argument("--min_called",    type=int,   default=5000)
    p.add_argument("--flag_low",      type=float, default=0.95,
                   help="Concordance below this flags same-patient pairs as mismatched")
    p.add_argument("--flag_high",     type=float, default=0.99,
                   help="Concordance above this flags different-patient pairs as matching")
    p.add_argument("--outdir",        default=".")
    return p.parse_args()


# ── VCF parsing ───────────────────────────────────────────────────────────────
def _get_info(info_str, tag):
    for field in info_str.split(";"):
        if field.startswith(f"{tag}="):
            return field[len(tag) + 1:]
    return None


def read_cellsnp_vcf(path, min_dp=20, max_oth_frac=0.10):
    rows = []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip().split('\t')
            if len(cols) < 8:
                continue
            chrom, pos, _, ref, alt = cols[:5]
            info = cols[7]
            ad_s  = _get_info(info, 'AD')
            dp_s  = _get_info(info, 'DP')
            oth_s = _get_info(info, 'OTH')
            if None in (ad_s, dp_s, oth_s):
                continue
            try:
                ad, dp, oth = float(ad_s), float(dp_s), float(oth_s)
            except ValueError:
                continue
            if dp < min_dp or oth > dp * max_oth_frac:
                continue
            if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                continue
            af = ad / dp
            if   af <= 0.10:            gt = 0
            elif 0.25 <= af <= 0.75:    gt = 1
            elif af >= 0.90:            gt = 2
            else:                       gt = None
            rows.append({'key': f"{chrom}:{pos}:{ref}:{alt}", 'gt': gt})
    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=['key', 'gt'])
    return df.dropna(subset=['gt']).reset_index(drop=True)


# ── Pairwise concordance ──────────────────────────────────────────────────────
def compare_pair(df_a, df_b):
    merged   = df_a.merge(df_b, on='key', suffixes=('_a', '_b'))
    n_shared = len(merged)
    called   = merged.dropna(subset=['gt_a', 'gt_b'])
    n_called = len(called)
    if n_called == 0:
        return n_shared, 0, None, None
    n_match = int((called['gt_a'] == called['gt_b']).sum())
    return n_shared, n_called, n_match, n_match / n_called


def interpret(concordance, n_called, min_called):
    if concordance is None or n_called < min_called:
        return "Too few sites — inconclusive"
    if concordance >= 0.99:   return "High — same donor"
    if concordance >= 0.95:   return "Intermediate — review"
    return "Low — likely different donors"


# ── Plot helpers ──────────────────────────────────────────────────────────────
def draw_flag_cell(ax, row_idx, col_idx, color, lw=4.0):
    """Make a flagged cell obvious: white halo + semi-transparent fill + solid border."""
    ax.add_patch(Rectangle(
        (col_idx, row_idx), 1, 1,
        fill=False, edgecolor='white', linewidth=lw + 2.5,
        zorder=5, clip_on=False
    ))
    ax.add_patch(Rectangle(
        (col_idx, row_idx), 1, 1,
        fill=True, facecolor=color, alpha=0.45,
        edgecolor='none', linewidth=0, zorder=6, clip_on=False
    ))
    ax.add_patch(Rectangle(
        (col_idx, row_idx), 1, 1,
        fill=False, edgecolor=color, linewidth=lw,
        zorder=7, clip_on=False
    ))


def hatch_inconclusive(ax, samples, nsites_df, min_called):
    """Overlay diagonal stripes on cells that have too few called sites."""
    for i, si in enumerate(samples):
        for j, sj in enumerate(samples):
            if i != j and int(nsites_df.loc[si, sj]) < min_called:
                # slight bleach so the colour still hints at direction
                ax.add_patch(Rectangle((j, i), 1, 1,
                    fill=True, facecolor='white', alpha=0.30,
                    edgecolor='none', zorder=3))
                # hatch lines
                ax.add_patch(Rectangle((j, i), 1, 1,
                    fill=False, hatch='//', edgecolor='#777777',
                    linewidth=0.3, zorder=4))


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    args   = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load samplesheet ──────────────────────────────────────────────────────
    ss = pd.read_excel(args.samplesheet)
    ss.columns = [c.strip().lower() for c in ss.columns]
    sid_col = next((c for c in ['sample_id','sample']           if c in ss.columns), None)
    pat_col = next((c for c in ['patient_id','patient','donor_id','donor'] if c in ss.columns), None)

    def get_patient(sid):
        if sid_col and pat_col:
            row = ss[ss[sid_col] == sid]
            if not row.empty:
                return str(row[pat_col].iloc[0])
        return 'unknown'

    if args.skip_genotyping:
        # ── Reload pre-computed results ───────────────────────────────────────
        tsv_path = outdir / "concordance_all_pairs.tsv"
        if not tsv_path.exists():
            sys.exit(f"ERROR: --skip_genotyping requires {tsv_path}\n"
                     f"       Run without --skip_genotyping first to generate it.")
        results = pd.read_csv(tsv_path, sep='\t').sort_values(
            ['same_patient', 'concordance'], ascending=[False, False])

        all_sids = sorted(set(results['sample_A']) | set(results['sample_B']))
        pat_groups = defaultdict(list)
        for sid in all_sids:
            pat_groups[get_patient(sid)].append(sid)
        sorted_samples = []
        for pat in sorted(pat_groups):
            sorted_samples.extend(sorted(pat_groups[pat]))
        n = len(sorted_samples)

        conc_matrix = pd.DataFrame(np.nan, index=sorted_samples, columns=sorted_samples)
        for i in range(n):
            conc_matrix.iloc[i, i] = 1.0
        for _, rr in results.iterrows():
            c = rr['concordance'] if pd.notna(rr['concordance']) else np.nan
            conc_matrix.loc[rr['sample_A'], rr['sample_B']] = c
            conc_matrix.loc[rr['sample_B'], rr['sample_A']] = c

        print(f"[compare_genotypes] --skip_genotyping: "
              f"loaded {len(results)} pairs from {tsv_path}")

    else:
        # ── Read VCFs ─────────────────────────────────────────────────────────
        if args.vcf_dir is None:
            sys.exit("ERROR: --vcf_dir is required unless --skip_genotyping is used")
        vcf_dir   = Path(args.vcf_dir)
        vcf_files = sorted(vcf_dir.glob("cellSNP.*.vcf.gz"))
        if not vcf_files:
            sys.exit(f"ERROR: no cellSNP.*.vcf.gz files in {vcf_dir}")

        print(f"[compare_genotypes] {len(vcf_files)} VCF files")
        sample_ids, vcf_data = [], {}
        for f in vcf_files:
            sid = re.sub(r'^cellSNP\.|\.vcf\.gz$', '', f.name)
            df  = read_cellsnp_vcf(f, min_dp=args.min_dp, max_oth_frac=args.max_oth_frac)
            print(f"  {sid}: {len(df):,} callable SNPs")
            sample_ids.append(sid)
            vcf_data[sid] = df

        n = len(sample_ids)

        # ── Sort samples by patient so same-patient pairs are adjacent ────────
        pat_groups = defaultdict(list)
        for sid in sample_ids:
            pat_groups[get_patient(sid)].append(sid)
        sorted_samples = []
        for pat in sorted(pat_groups):
            sorted_samples.extend(sorted(pat_groups[pat]))

        # ── Pairwise concordance matrix ───────────────────────────────────────
        conc_matrix = pd.DataFrame(np.nan, index=sorted_samples, columns=sorted_samples)
        for i in range(n):
            conc_matrix.iloc[i, i] = 1.0

        rows = []
        for i in range(n):
            for j in range(i + 1, n):
                a, b = sorted_samples[i], sorted_samples[j]
                n_shared, n_called, n_match, conc = compare_pair(vcf_data[a], vcf_data[b])
                conc_matrix.loc[a, b] = conc
                conc_matrix.loc[b, a] = conc
                pat_a, pat_b = get_patient(a), get_patient(b)
                same_pat     = (pat_a == pat_b) and (pat_a != 'unknown')
                rows.append({
                    'sample_A': a, 'sample_B': b,
                    'patient_A': pat_a, 'patient_B': pat_b,
                    'same_patient':   same_pat,
                    'shared_sites':   n_shared,
                    'called_sites':   n_called,
                    'matching_sites': n_match,
                    'concordance':    round(conc, 4) if conc is not None else None,
                    'interpretation': interpret(conc, n_called, args.min_called),
                })

        results = pd.DataFrame(rows).sort_values(
            ['same_patient', 'concordance'], ascending=[False, False])

        # ── Save tables ───────────────────────────────────────────────────────
        results.to_csv(outdir / "concordance_all_pairs.tsv",    sep='\t', index=False)
        results[results['same_patient']].to_csv(
            outdir / "concordance_same_patient.tsv", sep='\t', index=False)

    # ── Flagged pairs ─────────────────────────────────────────────────────────
    flag_intra = results[
        results['same_patient'] &
        results['concordance'].notna() &
        (results['called_sites'] >= args.min_called) &
        (results['concordance'] < args.flag_low)
    ].copy()
    flag_intra['flag'] = 'intra_patient_mismatch'

    flag_inter = results[
        ~results['same_patient'] &
        results['concordance'].notna() &
        (results['called_sites'] >= args.min_called) &
        (results['concordance'] >= args.flag_high)
    ].copy()
    flag_inter['flag'] = 'inter_patient_match'

    flags = pd.concat([flag_intra, flag_inter]).sort_values('concordance')

    if not args.skip_genotyping:
        flags.to_csv(outdir / "concordance_flags.tsv", sep='\t', index=False)

    print(f"\n  {len(results)} pairs  |  {len(flag_intra)} intra-mismatch  |  "
          f"{len(flag_inter)} inter-match flags")

    if len(flags):
        print("\n  ⚠ FLAGGED PAIRS:")
        for _, r in flags.iterrows():
            sym = "🔴" if r['flag'] == 'intra_patient_mismatch' else "🟠"
            print(f"  {sym} {r['sample_A']} × {r['sample_B']}  "
                  f"conc={r['concordance']}  n={r['called_sites']}  "
                  f"[{r['flag']}]")

    # ── Shared helpers ────────────────────────────────────────────────────────
    n_flags        = len(flag_intra) + len(flag_inter)

    # Dynamic colour scale — base vmin on actual data so values like 0.64 aren't
    # crushed to the bottom of a 0.90–1.00 scale
    _off = conc_matrix.values.astype(float).copy()
    np.fill_diagonal(_off, np.nan)
    _valid = _off[~np.isnan(_off)]
    _data_min = float(np.nanmin(_valid)) if len(_valid) else 0.70
    vmin = max(0.0, np.floor(_data_min * 20) / 20)   # round down to nearest 0.05
    vmax = 1.00
    flag_txt       = f"  —  ⚠ {n_flags} flagged pair(s)" if n_flags else "  —  no flags"
    idx            = {s: i for i, s in enumerate(sorted_samples)}

    # Called-sites matrix — needed for inconclusive hatching & n= annotation
    nsites_matrix = pd.DataFrame(0, index=sorted_samples, columns=sorted_samples)
    for s in sorted_samples:
        nsites_matrix.loc[s, s] = args.min_called + 1   # diagonal is never inconclusive
    for _, rr in results.iterrows():
        v = int(rr['called_sites']) if pd.notna(rr['called_sites']) else 0
        nsites_matrix.loc[rr['sample_A'], rr['sample_B']] = v
        nsites_matrix.loc[rr['sample_B'], rr['sample_A']] = v

    legend_handles = [
        mpatches.Patch(facecolor='#e31a1c',
                       label=f'Same patient, low concordance (<{args.flag_low})'),
        mpatches.Patch(facecolor='#ff7f00',
                       label=f'Diff patient, high concordance (≥{args.flag_high})'),
        mpatches.Patch(facecolor='white', hatch='//', edgecolor='#777777',
                       label=f'Too few sites (<{args.min_called:,} called)'),
    ]

    # Flagged samples for zoom page
    flag_samples_set = (set(flag_intra['sample_A']) | set(flag_intra['sample_B']) |
                        set(flag_inter['sample_A']) | set(flag_inter['sample_B']))
    flag_samples_ord = [s for s in sorted_samples if s in flag_samples_set]
    has_zoom = len(flag_samples_ord) >= 2
    nz       = len(flag_samples_ord)

    out_pdf = outdir / "genotype_concordance_heatmap.pdf"
    with PdfPages(out_pdf) as pdf:

        # ══════════════════════════════════════════════════════════════════════
        # PAGE 1 — full concordance heatmap
        # ══════════════════════════════════════════════════════════════════════
        cell_size = max(0.55, min(0.9, 20 / n))
        heat_side = n * cell_size
        cbar_in   = 2.8
        fig_w     = heat_side + cbar_in
        fig_h     = heat_side + 0.8

        fig1, ax_heat = plt.subplots(1, 1, figsize=(fig_w, fig_h))

        data = conc_matrix.loc[sorted_samples, sorted_samples].astype(float)
        mask = data.isna()

        # Annotation: show concordance + n_called only when matrix is small enough to read
        if n <= 35:
            ann1 = pd.DataFrame('', index=sorted_samples, columns=sorted_samples)
            for s in sorted_samples:
                ann1.loc[s, s] = '1.00'
            for _, rr in results.iterrows():
                a, b = rr['sample_A'], rr['sample_B']
                c = rr['concordance']
                nc = int(rr['called_sites']) if pd.notna(rr['called_sites']) else 0
                txt = (f"{c:.2f}\n{nc:,}" if pd.notna(c) else f"NA\n{nc:,}")
                ann1.loc[a, b] = txt
                ann1.loc[b, a] = txt
            annot_arg1 = ann1
            fmt_arg1   = ''
        else:
            annot_arg1 = False
            fmt_arg1   = '.2f'

        sns.heatmap(
            data, ax=ax_heat, mask=mask,
            cmap='RdYlGn', vmin=vmin, vmax=vmax,
            annot=annot_arg1, fmt=fmt_arg1,
            annot_kws={'size': max(5, 8 - n // 15), 'va': 'center'},
            linewidths=0.6, linecolor='#aaaaaa',
            square=True, xticklabels=True, yticklabels=True,
            cbar_kws={'label': 'Genotype concordance', 'shrink': 0.55, 'pad': 0.03},
        )

        fs = max(6, 10 - n // 10)
        ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=90, fontsize=fs)
        ax_heat.set_yticklabels(ax_heat.get_yticklabels(), rotation=0,  fontsize=fs)
        cbar1 = ax_heat.collections[0].colorbar
        cbar1.set_label('Genotype concordance', size=10)
        cbar1.ax.tick_params(labelsize=9)
        ax_heat.set_xlabel(""); ax_heat.set_ylabel("")

        # Patient dividers
        prev_pat, dividers = None, []
        for i, sid in enumerate(sorted_samples):
            pat = get_patient(sid)
            if pat != prev_pat and prev_pat is not None:
                dividers.append(i)
            prev_pat = pat
        for d in dividers:
            ax_heat.axhline(d, color='#444444', lw=1.5, zorder=4)
            ax_heat.axvline(d, color='#444444', lw=1.5, zorder=4)

        # Hatch cells with too few called sites (inconclusive)
        hatch_inconclusive(ax_heat, sorted_samples, nsites_matrix, args.min_called)

        # Flag borders
        for _, row in flag_intra.iterrows():
            i, j = idx[row['sample_A']], idx[row['sample_B']]
            for r, c in [(i, j), (j, i)]:
                draw_flag_cell(ax_heat, r, c, color='#e31a1c', lw=2.5)
        for _, row in flag_inter.iterrows():
            i, j = idx[row['sample_A']], idx[row['sample_B']]
            for r, c in [(i, j), (j, i)]:
                draw_flag_cell(ax_heat, r, c, color='#ff7f00', lw=2.5)

        ax_heat.legend(handles=legend_handles, loc='lower left',
                       bbox_to_anchor=(0.01, 0.01), ncol=1, fontsize=9, frameon=True,
                       framealpha=0.95, title='Flag types', title_fontsize=10)
        ax_heat.set_title(
            f"Cross-sample genotype concordance (cellsnp-lite pseudo-bulk){flag_txt}\n"
            f"Min DP ≥ {args.min_dp}  |  min called sites ≥ {args.min_called}",
            fontsize=9, pad=8)

        pdf.savefig(fig1, dpi=150, bbox_inches='tight')
        plt.close(fig1)

        # ══════════════════════════════════════════════════════════════════════
        # PAGE 2 — zoomed view of flagged samples only (only if flags exist)
        # ══════════════════════════════════════════════════════════════════════
        if has_zoom:
            zoom_cell  = max(0.9, min(2.5, 20 / nz))
            zoom_side  = nz * zoom_cell
            z_cbar_in  = 2.8
            z_fig_w    = zoom_side + z_cbar_in
            z_fig_h    = zoom_side + 1.0

            fig2, ax_z = plt.subplots(1, 1, figsize=(z_fig_w, z_fig_h))

            zoom_data = conc_matrix.loc[flag_samples_ord, flag_samples_ord].astype(float)
            zoom_mask = zoom_data.isna()
            zoom_idx  = {s: i for i, s in enumerate(flag_samples_ord)}

            # Annotation for zoom: concordance + n_called per cell
            ann2 = pd.DataFrame('', index=flag_samples_ord, columns=flag_samples_ord)
            for s in flag_samples_ord:
                ann2.loc[s, s] = '1.000'
            for _, rr in results.iterrows():
                a, b = rr['sample_A'], rr['sample_B']
                if a in zoom_idx and b in zoom_idx:
                    c  = rr['concordance']
                    nc = int(rr['called_sites']) if pd.notna(rr['called_sites']) else 0
                    txt = (f"{c:.3f}\nn={nc:,}" if pd.notna(c) else f"NA\nn={nc:,}")
                    ann2.loc[a, b] = txt
                    ann2.loc[b, a] = txt

            sns.heatmap(
                zoom_data, ax=ax_z, mask=zoom_mask,
                cmap='RdYlGn', vmin=vmin, vmax=vmax,
                annot=ann2, fmt='',
                annot_kws={'size': max(7, 12 - nz // 5), 'va': 'center'},
                linewidths=0.6, linecolor='#aaaaaa',
                square=True, xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'Genotype concordance', 'shrink': 0.55, 'pad': 0.03},
            )

            fs_z = max(9, 14 - nz // 4)
            ax_z.set_xticklabels(ax_z.get_xticklabels(), rotation=90, fontsize=fs_z)
            ax_z.set_yticklabels(ax_z.get_yticklabels(), rotation=0,  fontsize=fs_z)
            cbar2 = ax_z.collections[0].colorbar
            cbar2.set_label('Genotype concordance', size=11)
            cbar2.ax.tick_params(labelsize=10)
            ax_z.set_xlabel(""); ax_z.set_ylabel("")

            # Patient dividers
            prev_pat_z = None
            for i, sid in enumerate(flag_samples_ord):
                pat = get_patient(sid)
                if pat != prev_pat_z and prev_pat_z is not None:
                    ax_z.axhline(i, color='#444444', lw=2.0, zorder=4)
                    ax_z.axvline(i, color='#444444', lw=2.0, zorder=4)
                prev_pat_z = pat

            # Hatch inconclusive cells in zoom
            hatch_inconclusive(ax_z, flag_samples_ord, nsites_matrix, args.min_called)

            # Flag borders (thicker — cells are bigger here)
            for _, row in flag_intra.iterrows():
                if row['sample_A'] in zoom_idx and row['sample_B'] in zoom_idx:
                    i, j = zoom_idx[row['sample_A']], zoom_idx[row['sample_B']]
                    for r, c in [(i, j), (j, i)]:
                        draw_flag_cell(ax_z, r, c, color='#e31a1c', lw=4.0)
            for _, row in flag_inter.iterrows():
                if row['sample_A'] in zoom_idx and row['sample_B'] in zoom_idx:
                    i, j = zoom_idx[row['sample_A']], zoom_idx[row['sample_B']]
                    for r, c in [(i, j), (j, i)]:
                        draw_flag_cell(ax_z, r, c, color='#ff7f00', lw=4.0)

            ax_z.set_title(
                f"⚠  Flagged samples — zoomed view\n"
                f"{len(flag_intra)} intra-patient mismatch  |  "
                f"{len(flag_inter)} inter-patient match  |  "
                f"Min DP ≥ {args.min_dp}  |  min called sites ≥ {args.min_called}",
                fontsize=10, pad=10)

            pdf.savefig(fig2, dpi=150, bbox_inches='tight')
            plt.close(fig2)

    print(f"  Saved: {out_pdf}"
          + (f"  (2 pages: full heatmap + flagged zoom)" if has_zoom else ""))


if __name__ == "__main__":
    main()
