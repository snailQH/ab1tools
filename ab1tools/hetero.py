"""Heterozygous site detection and reporting.

Supports two modes:
    1. Frequency-based (original): uses base_freq dicts from BAM pileup
    2. SDVC signal-based (new): uses VariantCaller on chromatogram traces
"""

import csv
import os

BASES = "ACGT"

DEFAULT_MINOR_FREQ = 0.05   # minimum minor allele frequency to call hetero
DEFAULT_FLANKING = 20       # bases to show on each side of hetero site in PNG


def find_heterozygous_sites(base_freq, consensus_seq, min_minor_freq=DEFAULT_MINOR_FREQ):
    """Find positions with mixed/heterozygous signals.

    A position is heterozygous if the second-most-frequent base has
    frequency >= min_minor_freq.

    Returns list of dicts with keys:
        pos_0based, pos_1based, consensus_base, major_base, major_freq,
        minor_base, minor_freq, third_base, third_freq, freq_A, freq_C, freq_G, freq_T
    """
    sites = []

    for i, freq in enumerate(base_freq):
        sorted_bases = sorted(BASES, key=lambda b: freq.get(b, 0.0), reverse=True)
        major = sorted_bases[0]
        minor = sorted_bases[1]
        third = sorted_bases[2]
        major_freq = freq.get(major, 0.0)
        minor_freq = freq.get(minor, 0.0)
        third_freq = freq.get(third, 0.0)

        if minor_freq >= min_minor_freq:
            con_base = consensus_seq[i] if i < len(consensus_seq) else "N"
            sites.append({
                "pos_0based": i,
                "pos_1based": i + 1,
                "consensus_base": con_base,
                "major_base": major,
                "major_freq": major_freq,
                "minor_base": minor,
                "minor_freq": minor_freq,
                "third_base": third,
                "third_freq": third_freq,
                "freq_A": freq.get("A", 0.0),
                "freq_C": freq.get("C", 0.0),
                "freq_G": freq.get("G", 0.0),
                "freq_T": freq.get("T", 0.0),
            })

    return sites


def write_hetero_csv(sites, output_path, sample_name=""):
    """Write heterozygous sites report as CSV.

    Columns: sample, position, consensus, major_base, major_freq%,
             minor_base, minor_freq%, freq_A%, freq_C%, freq_G%, freq_T%
    """
    fieldnames = [
        "sample", "position", "consensus_base",
        "major_base", "major_freq_%",
        "minor_base", "minor_freq_%",
        "freq_A_%", "freq_C_%", "freq_G_%", "freq_T_%",
    ]

    write_header = not os.path.exists(output_path) or os.path.getsize(output_path) == 0

    with open(output_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()

        for site in sites:
            writer.writerow({
                "sample": sample_name,
                "position": site["pos_1based"],
                "consensus_base": site["consensus_base"],
                "major_base": site["major_base"],
                "major_freq_%": f"{site['major_freq'] * 100:.1f}",
                "minor_base": site["minor_base"],
                "minor_freq_%": f"{site['minor_freq'] * 100:.1f}",
                "freq_A_%": f"{site['freq_A'] * 100:.1f}",
                "freq_C_%": f"{site['freq_C'] * 100:.1f}",
                "freq_G_%": f"{site['freq_G'] * 100:.1f}",
                "freq_T_%": f"{site['freq_T'] * 100:.1f}",
            })

    return output_path


def get_hetero_windows(sites, n_bases, flanking=DEFAULT_FLANKING):
    """Compute display windows around heterozygous sites.

    Merges overlapping windows so nearby sites share one PNG panel.

    Returns list of (start_base, end_base) tuples (0-based, inclusive).
    """
    if not sites:
        return []

    # Build raw windows
    raw = []
    for site in sites:
        pos = site["pos_0based"]
        lo = max(0, pos - flanking)
        hi = min(n_bases - 1, pos + flanking)
        raw.append((lo, hi))

    # Merge overlapping windows
    raw.sort()
    merged = [raw[0]]
    for lo, hi in raw[1:]:
        prev_lo, prev_hi = merged[-1]
        if lo <= prev_hi + 1:
            merged[-1] = (prev_lo, max(prev_hi, hi))
        else:
            merged.append((lo, hi))

    return merged


def find_heterozygous_sites_sdvc(traces, peak_positions, base_calls,
                                  min_threshold=DEFAULT_MINOR_FREQ,
                                  sensitivity=0.5, prior=0.01):
    """Find heterozygous sites using the SDVC algorithm on trace data.

    This is the signal-level variant detection method that operates directly
    on 4-channel chromatogram traces, providing dynamic thresholding and
    multi-position consistency.

    Returns list of dicts with the same structure as find_heterozygous_sites()
    for backward compatibility, plus additional SDVC fields (confidence, snr).
    """
    from ab1tools.variant_caller import VariantCaller

    caller = VariantCaller(
        traces, peak_positions, base_calls,
        min_threshold=min_threshold,
        sensitivity=sensitivity,
        prior=prior,
    )
    all_results = caller.call_variants()

    sites = []
    for r in all_results:
        if not r["is_variant"]:
            continue
        sites.append({
            "pos_0based": r["pos_0based"],
            "pos_1based": r["pos_1based"],
            "consensus_base": r["consensus_base"],
            "major_base": r["major_base"],
            "major_freq": r["major_freq"],
            "minor_base": r["minor_base"],
            "minor_freq": r["minor_freq"],
            "third_base": r["third_base"],
            "third_freq": r["third_freq"],
            "freq_A": r["freq_A"],
            "freq_C": r["freq_C"],
            "freq_G": r["freq_G"],
            "freq_T": r["freq_T"],
            # SDVC-specific fields
            "confidence": r["confidence"],
            "snr": r["snr"],
            "threshold": r["threshold"],
        })

    return sites
