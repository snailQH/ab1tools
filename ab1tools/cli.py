"""CLI entry point for ab1tools."""

import argparse
import os
import sys
import glob

from ab1tools.signal import (
    extract_base_frequencies, read_consensus, generate_traces,
    read_sequence, sequence_to_base_frequencies,
    DEFAULT_SPACING, DEFAULT_SIGMA, DEFAULT_SCALE, DEFAULT_NOISE,
    DEFAULT_PHASING, DEFAULT_DECAY, DEFAULT_MIN_MAPQ,
)
from ab1tools.abif_writer import write_ab1
from ab1tools.abif_reader import read_ab1
from ab1tools.visualize import plot_chromatogram, plot_hetero_windows
from ab1tools.hetero import (
    find_heterozygous_sites, write_hetero_csv, get_hetero_windows,
    DEFAULT_MINOR_FREQ, DEFAULT_FLANKING,
)
from ab1tools.variant_caller import (
    VariantCaller, write_variants_csv, write_variants_vcf,
    subtract_control,
)


def convert_sample(bam_path, consensus_path, output_dir, sample_name=None,
                   spacing=DEFAULT_SPACING, sigma=DEFAULT_SIGMA,
                   scale=DEFAULT_SCALE, noise=DEFAULT_NOISE,
                   phasing=DEFAULT_PHASING, decay=DEFAULT_DECAY,
                   min_mapq=DEFAULT_MIN_MAPQ, use_consensus_calls=False,
                   no_plot=False, asymmetry=1.0, dye_scaling=None,
                   crosstalk=0.0, baseline_drift=0.0, smooth=0.0):
    """Convert a single sample (BAM + consensus) to AB1 (+ optional PNG)."""
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(bam_path))[0]

    os.makedirs(output_dir, exist_ok=True)

    ab1_path = os.path.join(output_dir, f"{sample_name}.ab1")
    png_path = os.path.join(output_dir, f"{sample_name}.png")

    # Read consensus
    print(f"  Reading consensus: {consensus_path}")
    consensus_seq, record = read_consensus(consensus_path)

    # Extract frequencies from BAM
    print(f"  Extracting base frequencies from: {bam_path}")
    base_freq = extract_base_frequencies(bam_path, consensus_seq, min_mapq=min_mapq)
    print(f"  Positions: {len(base_freq)}")

    # Generate traces
    print(f"  Generating traces (spacing={spacing}, sigma={sigma})")
    traces, base_calls, peak_positions, quality_scores = generate_traces(
        base_freq, spacing=spacing, sigma=sigma, scale=scale,
        noise=noise, phasing=phasing, decay=decay,
        asymmetry=asymmetry, dye_scaling=dye_scaling,
        crosstalk=crosstalk, baseline_drift=baseline_drift, smooth=smooth,
    )

    # Optionally override base calls with consensus sequence
    if use_consensus_calls:
        base_calls = consensus_seq[:len(base_calls)]
        print(f"  Using consensus base calls (overriding pileup calls)")

    # Write AB1
    write_ab1(ab1_path, traces, base_calls, peak_positions,
              sample_name=sample_name, quality_scores=quality_scores)
    print(f"  AB1 written: {ab1_path}")

    # Write PNGs (optional)
    png_paths = []
    if not no_plot:
        png_paths = plot_chromatogram(traces, base_calls, peak_positions, png_path,
                                     title=f"{sample_name} Chromatogram")
        for p in png_paths:
            print(f"  PNG written: {p}")

    return ab1_path, png_paths


def smart_sample(bam_path, consensus_path, output_dir, sample_name=None,
                 spacing=DEFAULT_SPACING, sigma=DEFAULT_SIGMA,
                 scale=DEFAULT_SCALE, noise=DEFAULT_NOISE,
                 phasing=DEFAULT_PHASING, decay=DEFAULT_DECAY,
                 min_mapq=DEFAULT_MIN_MAPQ, use_consensus_calls=False,
                 min_minor_freq=DEFAULT_MINOR_FREQ, flanking=DEFAULT_FLANKING,
                 csv_path=None, no_plot=False, asymmetry=1.0, dye_scaling=None,
                 crosstalk=0.0, baseline_drift=0.0, smooth=0.0):
    """Smart mode: AB1 + CSV report + PNGs only around heterozygous sites."""
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(bam_path))[0]

    os.makedirs(output_dir, exist_ok=True)

    ab1_path = os.path.join(output_dir, f"{sample_name}.ab1")

    # Read consensus
    print(f"  Reading consensus: {consensus_path}")
    consensus_seq, record = read_consensus(consensus_path)

    # Extract frequencies from BAM
    print(f"  Extracting base frequencies from: {bam_path}")
    base_freq = extract_base_frequencies(bam_path, consensus_seq, min_mapq=min_mapq)
    print(f"  Positions: {len(base_freq)}")

    # Generate traces
    print(f"  Generating traces (spacing={spacing}, sigma={sigma})")
    traces, base_calls, peak_positions, quality_scores = generate_traces(
        base_freq, spacing=spacing, sigma=sigma, scale=scale,
        noise=noise, phasing=phasing, decay=decay,
        asymmetry=asymmetry, dye_scaling=dye_scaling,
        crosstalk=crosstalk, baseline_drift=baseline_drift, smooth=smooth,
    )

    if use_consensus_calls:
        base_calls = consensus_seq[:len(base_calls)]

    # Write AB1
    write_ab1(ab1_path, traces, base_calls, peak_positions,
              sample_name=sample_name, quality_scores=quality_scores)
    print(f"  AB1 written: {ab1_path}")

    # Detect heterozygous sites
    hetero_sites = find_heterozygous_sites(base_freq, consensus_seq,
                                           min_minor_freq=min_minor_freq)
    print(f"  Heterozygous sites (minor >= {min_minor_freq*100:.0f}%): {len(hetero_sites)}")

    # Write to CSV
    if csv_path is None:
        csv_path = os.path.join(output_dir, "heterozygous_sites.csv")
    write_hetero_csv(hetero_sites, csv_path, sample_name=sample_name)

    # Generate PNGs only around hetero sites (optional)
    if hetero_sites and not no_plot:
        windows = get_hetero_windows(hetero_sites, len(base_freq), flanking=flanking)
        png_path = os.path.join(output_dir, f"{sample_name}_hetero.png")
        plot_hetero_windows(traces, base_calls, peak_positions,
                            hetero_sites, windows, png_path,
                            title=f"{sample_name}")
        print(f"  Hetero PNG written: {png_path}")
    if hetero_sites:
        for site in hetero_sites:
            print(f"    pos {site['pos_1based']}: {site['consensus_base']} -> "
                  f"{site['major_base']}={site['major_freq']*100:.1f}% / "
                  f"{site['minor_base']}={site['minor_freq']*100:.1f}%")
    else:
        print(f"  No heterozygous sites found — no PNG generated")

    return ab1_path, hetero_sites


def convert_sequence(input_path, output_dir, sample_name=None,
                     spacing=DEFAULT_SPACING, sigma=DEFAULT_SIGMA,
                     scale=DEFAULT_SCALE, noise=DEFAULT_NOISE,
                     phasing=DEFAULT_PHASING, decay=DEFAULT_DECAY,
                     no_plot=False, asymmetry=1.0, dye_scaling=None,
                     crosstalk=0.0, baseline_drift=0.0, smooth=0.0):
    """Convert a single FASTA/FASTQ file to AB1 + PNG (no BAM required).

    Generates clean single-peak traces from the sequence directly.
    """
    if sample_name is None:
        sample_name = os.path.splitext(os.path.basename(input_path))[0]

    os.makedirs(output_dir, exist_ok=True)

    ab1_path = os.path.join(output_dir, f"{sample_name}.ab1")
    png_path = os.path.join(output_dir, f"{sample_name}.png")

    # Read sequence (auto-detects FASTA or FASTQ)
    print(f"  Reading sequence: {input_path}")
    seq, record = read_sequence(input_path)
    print(f"  Sequence length: {len(seq)} bp")

    # Generate base frequencies from sequence (100% per base, no mixed signals)
    base_freq = sequence_to_base_frequencies(seq)

    # Generate traces
    print(f"  Generating traces (spacing={spacing}, sigma={sigma})")
    traces, base_calls, peak_positions, quality_scores = generate_traces(
        base_freq, spacing=spacing, sigma=sigma, scale=scale,
        noise=noise, phasing=phasing, decay=decay,
        asymmetry=asymmetry, dye_scaling=dye_scaling,
        crosstalk=crosstalk, baseline_drift=baseline_drift, smooth=smooth,
    )

    # Write AB1
    write_ab1(ab1_path, traces, base_calls, peak_positions,
              sample_name=sample_name, quality_scores=quality_scores)
    print(f"  AB1 written: {ab1_path}")

    # Write PNGs (optional)
    png_paths = []
    if not no_plot:
        png_paths = plot_chromatogram(traces, base_calls, peak_positions, png_path,
                                     title=f"{sample_name} Chromatogram")
        for p in png_paths:
            print(f"  PNG written: {p}")

    return ab1_path, png_paths


def plot_ab1_region(ab1_path, output_path, start=None, end=None,
                    bases_per_row=100, dpi=150):
    """Read an AB1 file and export PNG plots for a specific region.

    Args:
        ab1_path: path to input AB1 file
        output_path: path for output PNG file
        start: start base position (1-based, inclusive). None = beginning.
        end: end base position (1-based, inclusive). None = end of sequence.
        bases_per_row: number of bases to show per row in the PNG
        dpi: PNG resolution
    """
    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)

    n_bases = len(base_calls)
    if not n_bases:
        raise ValueError(f"No base calls found in {ab1_path}")

    # Convert 1-based inclusive to 0-based
    start_0 = (start - 1) if start else 0
    end_0 = (end - 1) if end else (n_bases - 1)

    # Clamp
    start_0 = max(0, min(start_0, n_bases - 1))
    end_0 = max(start_0, min(end_0, n_bases - 1))

    # Slice to the region of interest
    region_calls = base_calls[start_0:end_0 + 1]
    region_peaks = peak_positions[start_0:end_0 + 1]

    # Determine spacing from peak positions
    if len(peak_positions) > 1:
        spacing = peak_positions[1] - peak_positions[0]
    else:
        spacing = 10

    start_1 = start_0 + 1
    end_1 = end_0 + 1
    region_bases = end_1 - start_1 + 1

    print(f"  AB1 file: {ab1_path}")
    print(f"  Sample: {sample_name}")
    print(f"  Total bases: {n_bases}")
    print(f"  Plotting region: bases {start_1}–{end_1} ({region_bases} bp)")
    print(f"  Bases per row: {bases_per_row}")

    title = f"{sample_name} — Bases {start_1}–{end_1}"
    png_paths = plot_chromatogram(
        traces, base_calls, peak_positions, output_path,
        title=title, dpi=dpi, bases_per_row=bases_per_row,
        base_range=(start_0, end_0),
    )
    for p in png_paths:
        print(f"  PNG written: {p}")

    return png_paths


def view_ab1(ab1_path, as_json=False):
    """Display AB1 file metadata, trace statistics, and quality summary."""
    import json as json_mod
    import numpy as np

    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)
    n_bases = len(base_calls)
    trace_len = len(traces["A"]) if "A" in traces else 0
    spacing = (peak_positions[1] - peak_positions[0]) if len(peak_positions) > 1 else 0

    # Per-channel statistics
    channel_stats = {}
    for base in "ACGT":
        arr = traces[base].astype(float)
        channel_stats[base] = {
            "mean": round(float(np.mean(arr)), 1),
            "max": int(np.max(arr)),
            "min": int(np.min(arr)),
            "std": round(float(np.std(arr)), 1),
        }

    # Base composition
    comp = {b: base_calls.count(b) for b in "ACGT"}
    comp["N"] = n_bases - sum(comp.values())
    gc = (comp["G"] + comp["C"]) / max(n_bases, 1) * 100

    info = {
        "file": ab1_path,
        "sample_name": sample_name,
        "sequence_length": n_bases,
        "trace_length": trace_len,
        "peak_spacing": spacing,
        "base_composition": comp,
        "gc_content_%": round(gc, 1),
        "channel_stats": channel_stats,
    }

    if as_json:
        print(json_mod.dumps(info, indent=2))
    else:
        print(f"File:            {ab1_path}")
        print(f"Sample:          {sample_name}")
        print(f"Sequence length: {n_bases} bp")
        print(f"Trace length:    {trace_len} points")
        print(f"Peak spacing:    {spacing}")
        print(f"GC content:      {gc:.1f}%")
        print(f"Base composition: A={comp['A']} C={comp['C']} G={comp['G']} T={comp['T']}", end="")
        if comp["N"] > 0:
            print(f" N={comp['N']}")
        else:
            print()
        print(f"\nChannel statistics:")
        print(f"  {'Chan':>4s}  {'Mean':>7s}  {'Max':>6s}  {'Min':>6s}  {'Std':>7s}")
        for base in "ACGT":
            s = channel_stats[base]
            print(f"  {base:>4s}  {s['mean']:>7.1f}  {s['max']:>6d}  {s['min']:>6d}  {s['std']:>7.1f}")

    return info


def convert_ab1(ab1_path, output_format, output_path=None):
    """Convert AB1 file to FASTA, FASTQ, JSON, or CSV."""
    import json as json_mod
    import csv
    import numpy as np

    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)

    # Auto-generate output path
    if output_path is None:
        base = os.path.splitext(ab1_path)[0]
        ext_map = {"fasta": ".fasta", "fastq": ".fastq", "json": ".json", "csv": ".csv"}
        output_path = base + ext_map[output_format]

    if output_format == "fasta":
        with open(output_path, "w") as f:
            f.write(f">{sample_name}\n")
            # Wrap at 80 chars
            for i in range(0, len(base_calls), 80):
                f.write(base_calls[i:i+80] + "\n")

    elif output_format == "fastq":
        # Estimate quality from trace peak heights
        qualities = []
        for i, pos in enumerate(peak_positions):
            if i >= len(base_calls):
                break
            called_base = base_calls[i]
            if called_base in traces:
                arr = traces[called_base]
                lo = max(0, pos - 3)
                hi = min(len(arr), pos + 4)
                peak_val = float(np.max(arr[lo:hi]))
                # Simple quality estimation from peak height
                q = min(60, max(0, int(peak_val / 20)))
            else:
                q = 0
            qualities.append(q)

        with open(output_path, "w") as f:
            f.write(f"@{sample_name}\n")
            f.write(base_calls + "\n")
            f.write("+\n")
            f.write("".join(chr(q + 33) for q in qualities) + "\n")

    elif output_format == "json":
        data = {
            "sample_name": sample_name,
            "sequence": base_calls,
            "length": len(base_calls),
            "peak_positions": peak_positions[:len(base_calls)],
            "traces": {
                base: traces[base].tolist() for base in "ACGT" if base in traces
            },
        }
        with open(output_path, "w") as f:
            json_mod.dump(data, f, indent=2)

    elif output_format == "csv":
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "base", "peak_pos", "A", "C", "G", "T"])
            for i, pos in enumerate(peak_positions):
                if i >= len(base_calls):
                    break
                row = [i + 1, base_calls[i], pos]
                for base in "ACGT":
                    arr = traces[base]
                    lo = max(0, pos - 3)
                    hi = min(len(arr), pos + 4)
                    row.append(int(np.max(arr[lo:hi])))
                writer.writerow(row)

    print(f"  Converted: {ab1_path} → {output_path} ({output_format})")
    return output_path


def stats_ab1(ab1_path, output_path=None, output_format="text"):
    """Compute trace-level quality statistics from AB1 file."""
    import json as json_mod
    import csv
    import numpy as np

    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)
    n_bases = len(base_calls)

    # Per-position quality metrics
    positions = []
    for i, pos in enumerate(peak_positions):
        if i >= n_bases:
            break
        called = base_calls[i]
        signals = {}
        for base in "ACGT":
            arr = traces[base]
            lo = max(0, pos - 3)
            hi = min(len(arr), pos + 4)
            signals[base] = float(np.max(arr[lo:hi]))

        total = sum(signals.values())
        dominant = max(signals.values())
        noise_floor = min(signals.values())
        snr = dominant / max(noise_floor, 1.0)

        # Quality score from called base frequency
        freq = signals[called] / max(total, 1.0) if called in signals else 0
        if freq >= 1.0:
            q = 60
        elif freq <= 0:
            q = 0
        else:
            import math
            q = min(60, int(-10 * math.log10(max(1 - freq, 1e-10))))

        positions.append({
            "position": i + 1,
            "base": called,
            "quality": q,
            "snr": round(snr, 1),
            "peak_height": int(dominant),
            "signal_A": int(signals["A"]),
            "signal_C": int(signals["C"]),
            "signal_G": int(signals["G"]),
            "signal_T": int(signals["T"]),
        })

    # Summary statistics
    quals = [p["quality"] for p in positions]
    snrs = [p["snr"] for p in positions]
    peaks = [p["peak_height"] for p in positions]
    q20 = sum(1 for q in quals if q >= 20) / max(len(quals), 1) * 100
    q30 = sum(1 for q in quals if q >= 30) / max(len(quals), 1) * 100

    summary = {
        "sample": sample_name,
        "total_bases": n_bases,
        "mean_quality": round(float(np.mean(quals)), 1),
        "median_quality": int(np.median(quals)),
        "q20_%": round(q20, 1),
        "q30_%": round(q30, 1),
        "mean_snr": round(float(np.mean(snrs)), 1),
        "mean_peak_height": round(float(np.mean(peaks)), 1),
        "min_peak_height": int(np.min(peaks)),
        "max_peak_height": int(np.max(peaks)),
    }

    if output_format == "json":
        data = {"summary": summary, "per_position": positions}
        out = json_mod.dumps(data, indent=2)
    elif output_format == "csv":
        import io
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=positions[0].keys())
        writer.writeheader()
        for p in positions:
            writer.writerow(p)
        out = buf.getvalue()
    else:  # text
        lines = [
            f"Sample:           {summary['sample']}",
            f"Total bases:      {summary['total_bases']}",
            f"Mean quality:     Q{summary['mean_quality']}",
            f"Median quality:   Q{summary['median_quality']}",
            f"Q20 bases:        {summary['q20_%']}%",
            f"Q30 bases:        {summary['q30_%']}%",
            f"Mean SNR:         {summary['mean_snr']}",
            f"Mean peak height: {summary['mean_peak_height']}",
            f"Peak range:       {summary['min_peak_height']}-{summary['max_peak_height']}",
        ]
        out = "\n".join(lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(out + "\n")
        print(f"  Stats written: {output_path}")
    else:
        print(out)

    return summary


def extract_ab1_region(ab1_path, start=1, end=None, output_path=None):
    """Extract a sub-region from AB1 file into a new AB1 file.

    Like `samtools view -b file.bam chr1:100-200` but for AB1 files.
    Extracts base calls, traces, peak positions for the specified region
    and writes a new self-contained AB1 file.

    Args:
        ab1_path: path to input AB1 file
        start: start base position (1-based, inclusive)
        end: end base position (1-based, inclusive; None = end of sequence)
        output_path: path for output AB1 (auto-generated if None)
    """
    import numpy as np

    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)
    n = len(base_calls)

    if n == 0:
        print("  Warning: empty sequence, nothing to extract")
        return None

    # Convert 1-based to 0-based
    start_0 = max(0, start - 1)
    end_0 = min(n - 1, (end - 1) if end else (n - 1))

    if start_0 > end_0:
        print(f"  Error: start ({start}) > end ({end_0 + 1})")
        return None

    region_len = end_0 - start_0 + 1

    # Auto-generate output path
    if output_path is None:
        base, ext = os.path.splitext(ab1_path)
        output_path = f"{base}_{start}-{end_0 + 1}{ext}"

    # Extract base calls
    region_calls = base_calls[start_0:end_0 + 1]

    # Extract and rebase peak positions
    region_peaks_orig = peak_positions[start_0:end_0 + 1]
    if region_peaks_orig:
        trace_offset = region_peaks_orig[0]
        region_peaks = [p - trace_offset for p in region_peaks_orig]
    else:
        trace_offset = 0
        region_peaks = []

    # Determine spacing and trace bounds
    if len(region_peaks) > 1:
        spacing = region_peaks[1] - region_peaks[0]
    elif len(peak_positions) > 1:
        spacing = peak_positions[1] - peak_positions[0]
    else:
        spacing = 10

    # Extract trace region with padding on both sides
    trace_start = max(0, trace_offset - spacing)
    if region_peaks_orig:
        trace_end = region_peaks_orig[-1] + spacing
    else:
        trace_end = trace_offset + region_len * spacing

    # Adjust rebased peaks for the left padding
    left_pad = trace_offset - trace_start
    region_peaks = [p + left_pad for p in region_peaks]

    # Extract traces
    region_traces = {}
    for base in "ACGT":
        arr = traces[base]
        region_traces[base] = arr[trace_start:min(trace_end + 1, len(arr))]

    # Estimate quality scores from traces
    qualities = []
    for i, pos in enumerate(region_peaks):
        if i >= len(region_calls):
            break
        called = region_calls[i]
        if called in region_traces:
            arr = region_traces[called]
            lo = max(0, pos - 3)
            hi = min(len(arr), pos + 4)
            peak_val = float(np.max(arr[lo:hi]))
            total = 0
            for b in "ACGT":
                a = region_traces[b]
                total += float(np.max(a[lo:hi]))
            freq = peak_val / max(total, 1.0)
            if freq >= 1.0:
                q = 60
            elif freq <= 0:
                q = 0
            else:
                import math
                q = min(60, int(-10 * math.log10(max(1 - freq, 1e-10))))
        else:
            q = 0
        qualities.append(q)

    # Write new AB1
    region_name = f"{sample_name}_{start}-{end_0 + 1}"
    write_ab1(output_path, region_traces, region_calls, region_peaks,
              sample_name=region_name, quality_scores=qualities)

    print(f"  Input:    {ab1_path} ({n} bp)")
    print(f"  Region:   bases {start}–{end_0 + 1} ({region_len} bp)")
    print(f"  Traces:   {len(region_traces['A'])} points")
    print(f"  Output:   {output_path}")

    return {
        "input": ab1_path,
        "output": output_path,
        "original_length": n,
        "start": start,
        "end": end_0 + 1,
        "region_length": region_len,
        "trace_points": len(region_traces["A"]),
    }


def compare_ab1(ab1_path1, ab1_path2, output_path=None, output_format="text",
                plot=False, dpi=150):
    """Compare two AB1 files: signal correlation, base concordance, variant delta."""
    import json as json_mod
    import numpy as np

    t1, calls1, peaks1, name1 = read_ab1(ab1_path1)
    t2, calls2, peaks2, name2 = read_ab1(ab1_path2)

    n = min(len(calls1), len(calls2))

    # Base call concordance
    matches = sum(1 for i in range(n) if calls1[i] == calls2[i])
    concordance = matches / max(n, 1) * 100

    # Per-channel signal correlation (Pearson r)
    correlations = {}
    for base in "ACGT":
        a1 = t1[base].astype(np.float64)
        a2 = t2[base].astype(np.float64)
        min_len = min(len(a1), len(a2))
        if min_len > 0 and np.std(a1[:min_len]) > 0 and np.std(a2[:min_len]) > 0:
            r = float(np.corrcoef(a1[:min_len], a2[:min_len])[0, 1])
        else:
            r = 0.0
        correlations[base] = round(r, 4)

    mean_corr = round(float(np.mean(list(correlations.values()))), 4)

    # Mismatches detail
    mismatches = []
    for i in range(n):
        if calls1[i] != calls2[i]:
            mismatches.append({
                "position": i + 1,
                "file1_base": calls1[i],
                "file2_base": calls2[i],
            })

    result = {
        "file1": ab1_path1,
        "file2": ab1_path2,
        "sample1": name1,
        "sample2": name2,
        "compared_bases": n,
        "base_concordance_%": round(concordance, 2),
        "mismatches": len(mismatches),
        "channel_correlation": correlations,
        "mean_correlation": mean_corr,
        "mismatch_positions": mismatches[:50],  # cap at 50 for display
    }

    if output_format == "json":
        out = json_mod.dumps(result, indent=2)
    elif output_format == "csv":
        lines = ["metric,value"]
        lines.append(f"file1,{ab1_path1}")
        lines.append(f"file2,{ab1_path2}")
        lines.append(f"compared_bases,{n}")
        lines.append(f"base_concordance_%,{concordance:.2f}")
        lines.append(f"mismatches,{len(mismatches)}")
        lines.append(f"corr_A,{correlations['A']}")
        lines.append(f"corr_C,{correlations['C']}")
        lines.append(f"corr_G,{correlations['G']}")
        lines.append(f"corr_T,{correlations['T']}")
        lines.append(f"mean_correlation,{mean_corr}")
        out = "\n".join(lines)
    else:
        lines = [
            f"Comparing: {name1} vs {name2}",
            f"  File 1: {ab1_path1}",
            f"  File 2: {ab1_path2}",
            f"  Compared bases: {n}",
            f"",
            f"Base concordance: {concordance:.2f}% ({matches}/{n})",
            f"Mismatches:       {len(mismatches)}",
            f"",
            f"Signal correlation (Pearson r):",
            f"  A: {correlations['A']:.4f}",
            f"  C: {correlations['C']:.4f}",
            f"  G: {correlations['G']:.4f}",
            f"  T: {correlations['T']:.4f}",
            f"  Mean: {mean_corr:.4f}",
        ]
        if mismatches:
            lines.append(f"\nMismatch positions (first {min(20, len(mismatches))}):")
            for m in mismatches[:20]:
                lines.append(f"  pos {m['position']}: {m['file1_base']} → {m['file2_base']}")
        out = "\n".join(lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(out + "\n")
        print(f"  Comparison written: {output_path}")
    else:
        print(out)

    # Optional overlay plot
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        png_path = output_path.replace(".txt", ".png").replace(".csv", ".png") if output_path else "compare.png"
        if not png_path.endswith(".png"):
            png_path += ".png"

        fig, axes = plt.subplots(4, 1, figsize=(16, 8), sharex=True)
        colors = {"A": "#00CC00", "C": "#0000FF", "G": "#000000", "T": "#FF0000"}

        for idx, base in enumerate("ACGT"):
            ax = axes[idx]
            a1 = t1[base].astype(float)
            a2 = t2[base].astype(float)
            min_len = min(len(a1), len(a2))
            ax.plot(a1[:min_len], color=colors[base], alpha=0.7, linewidth=0.5, label=name1)
            ax.plot(a2[:min_len], color=colors[base], alpha=0.3, linewidth=0.5, linestyle="--", label=name2)
            ax.set_ylabel(base, fontsize=12, fontweight="bold")
            ax.legend(fontsize=8, loc="upper right")
            if idx == 0:
                ax.set_title(f"Overlay: {name1} vs {name2} (r={mean_corr:.4f})")

        axes[-1].set_xlabel("Trace position")
        plt.tight_layout()
        fig.savefig(png_path, dpi=dpi)
        plt.close(fig)
        print(f"  Overlay plot: {png_path}")

    return result


def trim_ab1(ab1_path, output_path=None, quality_threshold=20, min_length=50):
    """Trim low-quality regions from AB1 file using Mott's algorithm.

    Mott's algorithm (as implemented in phred, Ewing et al. 1998):
    1. Compute running score S_i = S_{i-1} + (quality_threshold - Q_i)
       starting from each end
    2. Find the region where the cumulative score is minimized on both ends
    3. The remaining region is the high-quality trimmed sequence
    """
    import numpy as np

    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)
    n = len(base_calls)

    if n == 0:
        print("  Warning: empty sequence, nothing to trim")
        return None

    # Estimate per-position quality from trace signals
    qualities = []
    for i, pos in enumerate(peak_positions):
        if i >= n:
            break
        called = base_calls[i]
        if called in traces:
            arr = traces[called]
            lo = max(0, pos - 3)
            hi = min(len(arr), pos + 4)
            peak_val = float(np.max(arr[lo:hi]))
            # All-channel total
            total = 0
            for b in "ACGT":
                a = traces[b]
                total += float(np.max(a[lo:hi]))
            freq = peak_val / max(total, 1.0)
            if freq >= 1.0:
                q = 60
            elif freq <= 0:
                q = 0
            else:
                import math
                q = min(60, int(-10 * math.log10(max(1 - freq, 1e-10))))
        else:
            q = 0
        qualities.append(q)

    # Mott's algorithm: find optimal trim region
    # Forward pass: find left trim point
    scores = [0.0]
    for q in qualities:
        scores.append(scores[-1] + (quality_threshold - q))

    # Find the region with maximum score drop (= best quality region)
    # The trim start is where cumulative score is maximum before the good region
    # The trim end is where cumulative score is minimum after the good region
    min_score = 0.0
    min_idx = 0
    max_drop = 0.0
    trim_start = 0
    trim_end = n - 1

    # Find best region using the max-subarray approach
    # Score = quality_threshold - Q_i; best region minimizes the sum
    # Equivalent to finding max contiguous sum of (Q_i - quality_threshold)
    adjusted = [q - quality_threshold for q in qualities]

    # Kadane's algorithm variant to find max-sum subarray
    best_sum = float("-inf")
    current_sum = 0
    current_start = 0

    for i, val in enumerate(adjusted):
        current_sum += val
        if current_sum > best_sum:
            best_sum = current_sum
            trim_start = current_start
            trim_end = i
        if current_sum < 0:
            current_sum = 0
            current_start = i + 1

    # Validate trim region
    trim_len = trim_end - trim_start + 1
    if trim_len < min_length:
        print(f"  Warning: trimmed region ({trim_len} bp) shorter than min_length ({min_length})")
        print(f"  Keeping full sequence")
        trim_start = 0
        trim_end = n - 1
        trim_len = n

    # Generate output path
    if output_path is None:
        base, ext = os.path.splitext(ab1_path)
        output_path = f"{base}_trimmed{ext}"

    # Extract trimmed data
    trimmed_calls = base_calls[trim_start:trim_end + 1]
    trimmed_peaks = peak_positions[trim_start:trim_end + 1]
    trimmed_quals = qualities[trim_start:trim_end + 1]

    # Rebase peak positions to start from 0
    if trimmed_peaks:
        offset = trimmed_peaks[0]
        trimmed_peaks = [p - offset for p in trimmed_peaks]

    # Determine spacing
    if len(trimmed_peaks) > 1:
        spacing = trimmed_peaks[1] - trimmed_peaks[0]
    else:
        spacing = 10

    # Extract trimmed traces
    if peak_positions and trim_start < len(peak_positions) and trim_end < len(peak_positions):
        trace_start = peak_positions[trim_start] - spacing
        trace_end = peak_positions[min(trim_end, len(peak_positions) - 1)] + spacing
        trace_start = max(0, trace_start)
    else:
        trace_start = 0
        trace_end = len(traces["A"])

    trimmed_traces = {}
    for base in "ACGT":
        arr = traces[base]
        trimmed_traces[base] = arr[trace_start:min(trace_end + 1, len(arr))]

    # Write trimmed AB1
    write_ab1(output_path, trimmed_traces, trimmed_calls, trimmed_peaks,
              sample_name=f"{sample_name}_trimmed", quality_scores=trimmed_quals)

    mean_q = float(np.mean(qualities[trim_start:trim_end + 1])) if trim_start <= trim_end else 0
    print(f"  Input:    {ab1_path} ({n} bp)")
    print(f"  Trimmed:  bases {trim_start + 1}–{trim_end + 1} ({trim_len} bp)")
    print(f"  Removed:  {trim_start} bp (left) + {n - trim_end - 1} bp (right)")
    print(f"  Mean Q:   {mean_q:.1f} (trimmed region)")
    print(f"  Output:   {output_path}")

    return {
        "input": ab1_path,
        "output": output_path,
        "original_length": n,
        "trim_start": trim_start + 1,
        "trim_end": trim_end + 1,
        "trimmed_length": trim_len,
        "removed_left": trim_start,
        "removed_right": n - trim_end - 1,
        "mean_quality": round(mean_q, 1),
    }


def call_ab1_variants(ab1_path, output_path=None, output_format="csv",
                      min_threshold=0.05, sensitivity=0.5, prior=0.01,
                      window_size=5, control_path=None,
                      plot=False, dpi=150):
    """Run SDVC variant calling on an AB1 file.

    Args:
        ab1_path: path to input AB1 file
        output_path: path for output file (CSV or VCF). None = auto-generate.
        output_format: "csv" or "vcf"
        min_threshold: minimum minor allele frequency (default 0.05)
        sensitivity: dynamic threshold sensitivity (default 0.5)
        prior: Bayesian prior for variant (default 0.01)
        window_size: multi-position consistency window (default 5)
        control_path: path to matched control AB1 file (Mode 2). None = Mode 1.
        plot: if True, generate annotated chromatogram PNG with variant highlights
        dpi: PNG resolution
    """
    traces, base_calls, peak_positions, sample_name = read_ab1(ab1_path)

    mode = "Control-free (Mode 1)"
    if control_path:
        control_traces, _, _, control_name = read_ab1(control_path)
        traces = subtract_control(traces, control_traces)
        mode = f"Control-enhanced (Mode 2, control={control_name})"

    print(f"  AB1 file: {ab1_path}")
    print(f"  Sample: {sample_name}")
    print(f"  Mode: {mode}")
    print(f"  Bases: {len(base_calls)}")

    caller = VariantCaller(
        traces, peak_positions, base_calls,
        min_threshold=min_threshold, sensitivity=sensitivity,
        prior=prior, window_size=window_size,
    )

    all_results = caller.call_variants()
    variants = [r for r in all_results if r["is_variant"]]

    print(f"  Variants detected: {len(variants)}")

    # Auto-generate output path if not provided
    if output_path is None:
        base = os.path.splitext(ab1_path)[0]
        ext = ".vcf" if output_format == "vcf" else ".variants.csv"
        output_path = f"{base}{ext}"

    # Write output
    if output_format == "vcf":
        write_variants_vcf(variants, output_path, sample_name=sample_name)
    else:
        write_variants_csv(variants, output_path, sample_name=sample_name)

    print(f"  Output: {output_path}")

    # Summary
    for v in variants:
        freq_str = f"{v['minor_freq']*100:.1f}%"
        conf_str = f"{v['confidence']:.3f}"
        print(f"    pos {v['pos_1based']}: {v['major_base']}→{v['minor_base']} "
              f"({freq_str}) conf={conf_str}")

    # Optional annotated plot
    if plot:
        from ab1tools.hetero import get_hetero_windows
        png_path = os.path.splitext(output_path)[0] + ".png"
        if variants:
            windows = get_hetero_windows(variants, len(base_calls))
            plot_hetero_windows(
                traces, base_calls, peak_positions,
                variants, windows, png_path,
                title=f"{sample_name} — SDVC Variants", dpi=dpi,
            )
            print(f"  Plot: {png_path}")

    return all_results, variants


def find_samples(analysis_dir):
    """Auto-detect barcode sample directories under an analysis directory."""
    samples = []
    for barcode_dir in sorted(glob.glob(os.path.join(analysis_dir, "barcode*"))):
        if not os.path.isdir(barcode_dir):
            continue
        barcode_name = os.path.basename(barcode_dir)

        bam_files = glob.glob(os.path.join(barcode_dir, "alignments", "*.aligned.sorted.bam"))
        if not bam_files:
            bam_files = glob.glob(os.path.join(barcode_dir, "alignments", "*.bam"))
        if not bam_files:
            print(f"  Warning: no BAM found in {barcode_dir}/alignments/, skipping")
            continue

        consensus_path = os.path.join(barcode_dir, "consensus", "consensus.fastq")
        if not os.path.exists(consensus_path):
            print(f"  Warning: no consensus.fastq in {barcode_dir}/consensus/, skipping")
            continue

        samples.append({
            "name": barcode_name,
            "bam": bam_files[0],
            "consensus": consensus_path,
        })

    return samples


def _add_signal_args(sub):
    """Add shared signal parameters to a subparser."""
    sub.add_argument("--spacing", type=int, default=DEFAULT_SPACING)
    sub.add_argument("--sigma", type=float, default=DEFAULT_SIGMA)
    sub.add_argument("--scale", type=int, default=DEFAULT_SCALE)
    sub.add_argument("--noise", type=int, default=DEFAULT_NOISE)
    sub.add_argument("--phasing", type=float, default=DEFAULT_PHASING)
    sub.add_argument("--decay", type=float, default=DEFAULT_DECAY,
                     help="Signal decay factor (0=uniform, 0.5=mild, 2.0=strong)")
    sub.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ)
    sub.add_argument("--use-consensus-calls", action="store_true",
                     help="Use consensus sequence for base calls instead of pileup majority")
    # Phase 1 signal realism parameters
    sub.add_argument("--asymmetry", type=float, default=1.0,
                     help="Peak asymmetry ratio (1.0=symmetric, 1.3=realistic rightward tail)")
    sub.add_argument("--dye-scaling", type=str, default=None,
                     help="Dye-specific scaling G,A,T,C (e.g., '0.80,1.00,0.95,0.85')")
    sub.add_argument("--crosstalk", type=float, default=0.0,
                     help="Cross-channel spectral crosstalk fraction (0=none, 0.03=realistic)")
    sub.add_argument("--baseline-drift", type=float, default=0.0,
                     help="Baseline drift amplitude (0=none, 15=mild)")
    sub.add_argument("--smooth", type=float, default=0.0,
                     help="Post-generation Gaussian smoothing sigma (0=none, 1.0=mild)")


def main():
    parser = argparse.ArgumentParser(
        prog="ab1tools",
        description="Comprehensive toolkit for AB1 chromatogram files: "
                    "generate, visualize, analyze, and convert — the samtools for AB1",
        epilog="Author: Qinghui Li (leeqhui@gmail.com) | "
               "https://github.com/snailQH/ab1tools",
    )
    subparsers = parser.add_subparsers(dest="command")

    # Single sample mode
    single = subparsers.add_parser("single", help="Convert a single sample")
    single.add_argument("--bam", required=True, help="Input BAM file")
    single.add_argument("--consensus", required=True, help="Input consensus FASTQ")
    single.add_argument("--output-dir", "-o", default=".", help="Output directory")
    single.add_argument("--name", help="Sample name (default: from BAM filename)")
    single.add_argument("--no-plot", action="store_true", default=False,
                        help="Skip PNG plot generation (AB1 only)")
    _add_signal_args(single)

    # Batch mode
    batch = subparsers.add_parser("batch", help="Batch convert all barcode samples")
    batch.add_argument("analysis_dir", help="Analysis directory containing barcodeXX/ subdirs")
    batch.add_argument("--output-dir", "-o", default="output", help="Output directory")
    batch.add_argument("--no-plot", action="store_true", default=False,
                       help="Skip PNG plot generation (AB1 only)")
    _add_signal_args(batch)

    # From-seq mode — FASTA/FASTQ only, no BAM needed
    fromseq = subparsers.add_parser("from-seq",
                                     help="Generate AB1 + PNG from a single FASTA or FASTQ file (no BAM required)")
    fromseq.add_argument("input", help="Input FASTA or FASTQ file")
    fromseq.add_argument("--output-dir", "-o", default=".", help="Output directory")
    fromseq.add_argument("--name", help="Sample name (default: from input filename)")
    fromseq.add_argument("--spacing", type=int, default=DEFAULT_SPACING)
    fromseq.add_argument("--sigma", type=float, default=DEFAULT_SIGMA)
    fromseq.add_argument("--scale", type=int, default=DEFAULT_SCALE)
    fromseq.add_argument("--noise", type=int, default=DEFAULT_NOISE)
    fromseq.add_argument("--phasing", type=float, default=DEFAULT_PHASING)
    fromseq.add_argument("--decay", type=float, default=DEFAULT_DECAY,
                         help="Signal decay factor (0=uniform, 0.5=mild, 2.0=strong)")
    fromseq.add_argument("--no-plot", action="store_true", default=False,
                         help="Skip PNG plot generation (AB1 only)")
    # Phase 1 signal realism parameters for from-seq
    fromseq.add_argument("--asymmetry", type=float, default=1.0,
                         help="Peak asymmetry ratio (1.0=symmetric, 1.3=realistic)")
    fromseq.add_argument("--dye-scaling", type=str, default=None,
                         help="Dye-specific scaling G,A,T,C (e.g., '0.80,1.00,0.95,0.85')")
    fromseq.add_argument("--crosstalk", type=float, default=0.0,
                         help="Cross-channel crosstalk (0=none, 0.03=realistic)")
    fromseq.add_argument("--baseline-drift", type=float, default=0.0,
                         help="Baseline drift amplitude (0=none, 15=mild)")
    fromseq.add_argument("--smooth", type=float, default=0.0,
                         help="Gaussian smoothing sigma (0=none, 1.0=mild)")

    # Plot-ab1 mode — export PNG from existing AB1 file with custom region
    plotab1 = subparsers.add_parser("plot-ab1",
                                    help="Export PNG chromatogram from an AB1 file with custom region and bin size")
    plotab1.add_argument("input", help="Input AB1 file")
    plotab1.add_argument("--output", "-o", help="Output PNG path (default: <input>.png)")
    plotab1.add_argument("--start", type=int, default=None,
                         help="Start base position (1-based, inclusive; default: beginning)")
    plotab1.add_argument("--end", type=int, default=None,
                         help="End base position (1-based, inclusive; default: end of sequence)")
    plotab1.add_argument("--bases-per-row", type=int, default=100,
                         help="Number of bases per row in the PNG (default: 100)")
    plotab1.add_argument("--dpi", type=int, default=150,
                         help="PNG resolution (default: 150)")

    # Call mode — SDVC variant calling on AB1 file
    call = subparsers.add_parser("call",
                                  help="Detect variants from AB1 file using SDVC algorithm")
    call.add_argument("input", help="Input AB1 file")
    call.add_argument("--output", "-o", help="Output file path (default: <input>.variants.csv)")
    call.add_argument("--format", choices=["csv", "vcf"], default="csv",
                      help="Output format: csv or vcf (default: csv)")
    call.add_argument("--threshold", type=float, default=0.05,
                      help="Minimum minor allele frequency threshold (default: 0.05)")
    call.add_argument("--sensitivity", type=float, default=0.5,
                      help="Dynamic threshold sensitivity: lower=more sensitive (default: 0.5)")
    call.add_argument("--prior", type=float, default=0.01,
                      help="Bayesian prior for variant probability (default: 0.01)")
    call.add_argument("--window-size", type=int, default=5,
                      help="Multi-position consistency window size (default: 5)")
    call.add_argument("--control", default=None,
                      help="Control AB1 file for enhanced detection (Mode 2, LOD <1%%)")
    call.add_argument("--plot", action="store_true", default=False,
                      help="Generate annotated chromatogram PNG with variant highlights")
    call.add_argument("--dpi", type=int, default=150,
                      help="PNG resolution for --plot (default: 150)")

    # View mode — inspect AB1 metadata and trace summary
    view = subparsers.add_parser("view",
                                  help="Display AB1 file metadata, trace stats, and quality summary")
    view.add_argument("input", help="Input AB1 file")
    view.add_argument("--json", action="store_true", default=False,
                      help="Output as JSON instead of human-readable text")

    # Convert mode — AB1 to other formats
    convert = subparsers.add_parser("convert",
                                     help="Convert AB1 to FASTA, FASTQ, JSON, or CSV")
    convert.add_argument("input", help="Input AB1 file")
    convert.add_argument("--format", "-f", required=True,
                         choices=["fasta", "fastq", "json", "csv"],
                         help="Output format: fasta, fastq, json, or csv")
    convert.add_argument("--output", "-o", help="Output file path (default: auto-generated)")

    # Stats mode — trace-level quality statistics
    stats = subparsers.add_parser("stats",
                                   help="Compute trace-level quality statistics from AB1 file")
    stats.add_argument("input", help="Input AB1 file")
    stats.add_argument("--output", "-o", help="Output file path (default: stdout)")
    stats.add_argument("--format", choices=["text", "csv", "json"], default="text",
                       help="Output format (default: text)")

    # Compare mode — compare two AB1 files
    compare = subparsers.add_parser("compare",
                                     help="Compare two AB1 files: signal correlation, base concordance, variant delta")
    compare.add_argument("input1", help="First AB1 file (e.g., sample)")
    compare.add_argument("input2", help="Second AB1 file (e.g., control or replicate)")
    compare.add_argument("--output", "-o", help="Output report file (default: stdout)")
    compare.add_argument("--format", choices=["text", "csv", "json"], default="text",
                         help="Output format (default: text)")
    compare.add_argument("--plot", action="store_true", default=False,
                         help="Generate overlay chromatogram PNG")
    compare.add_argument("--dpi", type=int, default=150,
                         help="PNG resolution for --plot (default: 150)")

    # Extract mode — extract sub-region from AB1 into new AB1
    extract = subparsers.add_parser("extract",
                                     help="Extract a sub-region from AB1 file into a new AB1 file")
    extract.add_argument("input", help="Input AB1 file")
    extract.add_argument("--start", "-s", type=int, default=1,
                         help="Start base position (1-based, inclusive; default: 1)")
    extract.add_argument("--end", "-e", type=int, default=None,
                         help="End base position (1-based, inclusive; default: end of sequence)")
    extract.add_argument("--output", "-o", help="Output AB1 file (default: <input>_<start>-<end>.ab1)")

    # Trim mode — quality-based trace trimming
    trim = subparsers.add_parser("trim",
                                  help="Trim low-quality regions from AB1 file (Mott's algorithm)")
    trim.add_argument("input", help="Input AB1 file")
    trim.add_argument("--output", "-o", help="Output trimmed AB1 file (default: <input>_trimmed.ab1)")
    trim.add_argument("--quality-threshold", type=int, default=20,
                      help="Minimum quality score threshold (default: Q20)")
    trim.add_argument("--min-length", type=int, default=50,
                      help="Minimum trimmed sequence length (default: 50)")

    # Smart mode — AB1 + hetero CSV + PNGs only at hetero sites
    smart = subparsers.add_parser("smart",
                                  help="Smart mode: AB1 + hetero report CSV + PNGs only at variant sites")
    smart.add_argument("analysis_dir", help="Analysis directory containing barcodeXX/ subdirs")
    smart.add_argument("--output-dir", "-o", default="output", help="Output directory")
    smart.add_argument("--no-plot", action="store_true", default=False,
                       help="Skip PNG plot generation (AB1 + CSV only)")
    smart.add_argument("--min-minor-freq", type=float, default=DEFAULT_MINOR_FREQ,
                       help=f"Minimum minor allele frequency to call heterozygous (default: {DEFAULT_MINOR_FREQ})")
    smart.add_argument("--flanking", type=int, default=DEFAULT_FLANKING,
                       help=f"Bases to show on each side of hetero site (default: {DEFAULT_FLANKING})")
    _add_signal_args(smart)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "view":
        view_ab1(args.input, as_json=args.json)
        return

    if args.command == "convert":
        convert_ab1(args.input, args.format, output_path=args.output)
        return

    if args.command == "stats":
        stats_ab1(args.input, output_path=args.output, output_format=args.format)
        return

    if args.command == "extract":
        extract_ab1_region(args.input, start=args.start, end=args.end,
                           output_path=args.output)
        return

    if args.command == "compare":
        compare_ab1(args.input1, args.input2,
                    output_path=args.output, output_format=args.format,
                    plot=args.plot, dpi=args.dpi)
        return

    if args.command == "trim":
        trim_ab1(args.input, output_path=args.output,
                 quality_threshold=args.quality_threshold,
                 min_length=args.min_length)
        return

    if args.command == "call":
        print(f"Running SDVC variant calling...")
        call_ab1_variants(
            args.input, output_path=args.output,
            output_format=args.format,
            min_threshold=args.threshold,
            sensitivity=args.sensitivity,
            prior=args.prior,
            window_size=args.window_size,
            control_path=args.control,
            plot=args.plot, dpi=args.dpi,
        )
        print("Done!")
        return

    if args.command == "from-seq":
        print(f"Converting sequence file to AB1...")
        # Parse dye-scaling for from-seq mode
        fs_dye = None
        if hasattr(args, 'dye_scaling') and args.dye_scaling:
            vals = [float(v) for v in args.dye_scaling.split(",")]
            if len(vals) == 4:
                fs_dye = {"G": vals[0], "A": vals[1], "T": vals[2], "C": vals[3]}
        convert_sequence(
            args.input, args.output_dir,
            sample_name=args.name,
            spacing=args.spacing, sigma=args.sigma, scale=args.scale,
            noise=args.noise, phasing=args.phasing, decay=args.decay,
            no_plot=args.no_plot,
            asymmetry=args.asymmetry, dye_scaling=fs_dye,
            crosstalk=args.crosstalk, baseline_drift=args.baseline_drift,
            smooth=args.smooth,
        )
        print("Done!")
        return

    if args.command == "plot-ab1":
        output = args.output
        if output is None:
            base = os.path.splitext(args.input)[0]
            suffix = ""
            if args.start or args.end:
                s = args.start or 1
                e = args.end or "end"
                suffix = f"_{s}-{e}"
            output = f"{base}{suffix}.png"
        print(f"Exporting chromatogram PNG from AB1...")
        plot_ab1_region(
            args.input, output,
            start=args.start, end=args.end,
            bases_per_row=args.bases_per_row, dpi=args.dpi,
        )
        print("Done!")
        return

    # Parse dye-scaling string if provided
    dye_scaling = None
    if hasattr(args, 'dye_scaling') and args.dye_scaling:
        vals = [float(v) for v in args.dye_scaling.split(",")]
        if len(vals) == 4:
            dye_scaling = {"G": vals[0], "A": vals[1], "T": vals[2], "C": vals[3]}

    signal_kwargs = dict(
        spacing=args.spacing, sigma=args.sigma, scale=args.scale,
        noise=args.noise, phasing=args.phasing, decay=args.decay,
        min_mapq=args.min_mapq, use_consensus_calls=args.use_consensus_calls,
        asymmetry=args.asymmetry,
        dye_scaling=dye_scaling,
        crosstalk=args.crosstalk,
        baseline_drift=args.baseline_drift,
        smooth=args.smooth,
    )

    if args.command == "single":
        print(f"Converting single sample...")
        convert_sample(
            args.bam, args.consensus, args.output_dir,
            sample_name=args.name, no_plot=args.no_plot, **signal_kwargs,
        )
        print("Done!")

    elif args.command == "batch":
        print(f"Scanning {args.analysis_dir} for samples...")
        samples = find_samples(args.analysis_dir)
        if not samples:
            print("No samples found!")
            sys.exit(1)

        print(f"Found {len(samples)} samples")
        for sample in samples:
            print(f"\n[{sample['name']}]")
            convert_sample(
                sample["bam"], sample["consensus"], args.output_dir,
                sample_name=sample["name"], no_plot=args.no_plot, **signal_kwargs,
            )

        print(f"\nAll done! Output in {args.output_dir}/")

    elif args.command == "smart":
        print(f"Smart mode: scanning {args.analysis_dir} for samples...")
        samples = find_samples(args.analysis_dir)
        if not samples:
            print("No samples found!")
            sys.exit(1)

        os.makedirs(args.output_dir, exist_ok=True)
        csv_path = os.path.join(args.output_dir, "heterozygous_sites.csv")

        # Remove old CSV so we start fresh
        if os.path.exists(csv_path):
            os.remove(csv_path)

        print(f"Found {len(samples)} samples")
        total_hetero = 0
        for sample in samples:
            print(f"\n[{sample['name']}]")
            _, sites = smart_sample(
                sample["bam"], sample["consensus"], args.output_dir,
                sample_name=sample["name"],
                min_minor_freq=args.min_minor_freq,
                flanking=args.flanking,
                csv_path=csv_path,
                no_plot=args.no_plot,
                **signal_kwargs,
            )
            total_hetero += len(sites)

        print(f"\nAll done! Output in {args.output_dir}/")
        print(f"  CSV report: {csv_path}")
        print(f"  Total heterozygous sites: {total_hetero}")


if __name__ == "__main__":
    main()
