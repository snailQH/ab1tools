"""Core signal generation: BAM pileup → 4-channel Sanger-like chromatogram traces.

Signal realism features (Phase 1):
  - Asymmetric peak modeling (EMG) [Grushka, Anal. Chem. 1972]
  - Dye-specific amplitude scaling [Ju et al., PNAS 1995]
  - Cross-channel spectral crosstalk [Huang et al., Anal. Chem. 1992]
  - Baseline drift (low-frequency wandering)
  - Local Gaussian smoothing [SciPy; Virtanen et al., Nat. Methods 2020]
  - Realistic S/N% computation [Giddings et al., Genome Res. 1998]
"""

import math
import numpy as np
import pysam
from Bio import SeqIO
from collections import Counter

BASES = "ACGT"

# Default signal parameters
DEFAULT_SPACING = 10    # points between peaks
DEFAULT_SIGMA = 2.0     # Gaussian peak width
DEFAULT_SCALE = 1024    # peak amplitude
DEFAULT_NOISE = 0       # random noise level (0 = clean)
DEFAULT_PHASING = 0.0   # tailing/drag effect (0 = none)
DEFAULT_DECAY = 0.0     # signal decay factor (0 = uniform)
DEFAULT_MIN_MAPQ = 20   # minimum mapping quality

# Phase 1 defaults
DEFAULT_ASYMMETRY = 1.0     # peak asymmetry ratio (1.0 = symmetric, 1.3 = realistic rightward tail)
DEFAULT_DYE_SCALING = None  # per-channel scale factors (None = uniform; realistic: G=0.80, A=1.00, T=0.95, C=0.85)
DEFAULT_CROSSTALK = 0.0     # spectral crosstalk fraction (0 = none, 0.03 = realistic)
DEFAULT_BASELINE_DRIFT = 0.0  # baseline drift amplitude (0 = none, 15 = mild)
DEFAULT_SMOOTH = 0.0        # post-generation Gaussian smoothing sigma (0 = none, 1.0 = mild)


def extract_base_frequencies(bam_path, consensus_seq, min_mapq=DEFAULT_MIN_MAPQ):
    """Extract per-position base frequencies from BAM pileup.

    Returns list of dicts like [{"A": 0.7, "C": 0.0, "G": 0.3, "T": 0.0}, ...]
    Falls back to consensus base when no coverage at a position.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Get reference name (first reference in BAM)
    ref_names = bam.references
    if not ref_names:
        raise ValueError(f"BAM file {bam_path} has no references")
    ref_name = ref_names[0]
    ref_len = bam.get_reference_length(ref_name)

    # Initialize with consensus fallback for all positions
    n = min(len(consensus_seq), ref_len)
    base_freq = []
    for i in range(n):
        freq = {b: 0.0 for b in BASES}
        if i < len(consensus_seq) and consensus_seq[i].upper() in BASES:
            freq[consensus_seq[i].upper()] = 1.0
        else:
            # Ambiguous base: distribute equally across ACGT
            for b in BASES:
                freq[b] = 0.25
        base_freq.append(freq)

    # Overwrite with actual pileup data where available
    for col in bam.pileup(ref_name, 0, n, stepper="all", truncate=True,
                          min_mapping_quality=min_mapq):
        pos = col.reference_pos
        if pos >= n:
            break

        counts = Counter()
        for read in col.pileups:
            if read.is_del or read.is_refskip:
                continue
            if read.query_position is None:
                continue
            base = read.alignment.query_sequence[read.query_position]
            if base in BASES:
                counts[base] += 1

        total = sum(counts.values())
        if total > 0:
            base_freq[pos] = {b: counts.get(b, 0) / total for b in BASES}

    bam.close()
    return base_freq


def read_consensus(file_path):
    """Read consensus sequence from FASTA or FASTQ file (auto-detected)."""
    fmt = detect_sequence_format(file_path)
    record = SeqIO.read(file_path, fmt)
    return str(record.seq), record


def detect_sequence_format(file_path):
    """Detect whether a file is FASTA or FASTQ based on the first character."""
    with open(file_path) as f:
        first_char = f.read(1)
    if first_char == ">":
        return "fasta"
    elif first_char == "@":
        return "fastq"
    else:
        raise ValueError(f"Cannot detect format of {file_path}: "
                         f"expected '>' (FASTA) or '@' (FASTQ), got '{first_char}'")


def read_sequence(file_path):
    """Read a sequence from a FASTA or FASTQ file (auto-detected).

    Returns:
        seq: str of the sequence
        record: Bio.SeqRecord
    """
    fmt = detect_sequence_format(file_path)
    record = SeqIO.read(file_path, fmt)
    return str(record.seq), record


def sequence_to_base_frequencies(seq):
    """Convert a plain sequence string to base frequencies (100% for each base).

    This is used when no BAM file is available — each position gets 100%
    frequency for its base, producing clean single-peak traces.

    Returns list of dicts like [{"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0}, ...]
    """
    base_freq = []
    for ch in seq.upper():
        freq = {b: 0.0 for b in BASES}
        if ch in BASES:
            freq[ch] = 1.0
        else:
            # Ambiguous base: distribute equally
            for b in BASES:
                freq[b] = 0.25
        base_freq.append(freq)
    return base_freq


def _gaussian(x, mu, sigma):
    """Symmetric Gaussian peak function."""
    return np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


def _asymmetric_gaussian(x, mu, sigma, asymmetry=1.0):
    """Asymmetric Gaussian (exponentially modified Gaussian approximation).

    Real Sanger peaks have rightward tailing due to band broadening in
    capillary electrophoresis. The asymmetry parameter controls the ratio
    of right-side to left-side sigma.

    Args:
        x: array of trace positions
        mu: peak center
        sigma: left-side standard deviation
        asymmetry: ratio of sigma_right / sigma_left (1.0 = symmetric,
                   1.3 = realistic rightward tail; Grushka, 1972)
    """
    sigma_left = sigma
    sigma_right = sigma * asymmetry
    result = np.where(
        x <= mu,
        np.exp(-(x - mu) ** 2 / (2 * sigma_left ** 2)),
        np.exp(-(x - mu) ** 2 / (2 * sigma_right ** 2)),
    )
    return result


def compute_snr(traces):
    """Compute per-channel signal-to-noise ratios from trace data.

    Returns dict with S/N% values for G, A, T, C channels in the order
    expected by the AB1 S/N%1 tag.

    SNR = mean(peak_heights) / std(inter-peak_baseline)
    """
    snr = {}
    for base in BASES:
        arr = traces[base].astype(np.float64)
        if len(arr) == 0:
            snr[base] = 0
            continue
        # Signal: 90th percentile of positive values (peak region)
        positive = arr[arr > 0]
        if len(positive) == 0:
            snr[base] = 0
            continue
        signal = float(np.percentile(positive, 90))
        # Noise: std of values below 20th percentile (baseline region)
        threshold = float(np.percentile(arr, 20))
        baseline = arr[arr <= max(threshold, 1)]
        noise_std = float(np.std(baseline)) if len(baseline) > 10 else 1.0
        snr[base] = int(min(999, signal / max(noise_std, 1.0)))
    return snr


def generate_traces(base_freq, spacing=DEFAULT_SPACING, sigma=DEFAULT_SIGMA,
                    scale=DEFAULT_SCALE, noise=DEFAULT_NOISE,
                    phasing=DEFAULT_PHASING, decay=DEFAULT_DECAY,
                    asymmetry=DEFAULT_ASYMMETRY,
                    dye_scaling=DEFAULT_DYE_SCALING,
                    crosstalk=DEFAULT_CROSSTALK,
                    baseline_drift=DEFAULT_BASELINE_DRIFT,
                    smooth=DEFAULT_SMOOTH):
    """Generate 4-channel chromatogram traces from base frequencies.

    Args:
        base_freq: list of dicts with ACGT frequencies per position
        spacing: points between peaks (default 10)
        sigma: Gaussian peak width (default 2.0)
        scale: max peak amplitude (default 1024)
        noise: noise level (0 = clean, 5 = realistic)
        phasing: tailing effect (0 = none, 0.1 = realistic)
        decay: signal decay (0 = uniform, 0.5 = mild, 2.0 = strong)
        asymmetry: peak asymmetry ratio (1.0 = symmetric, 1.3 = realistic
                   rightward tail; Grushka 1972)
        dye_scaling: per-channel scale factors dict, e.g.
                     {"G": 0.80, "A": 1.00, "T": 0.95, "C": 0.85}
                     (None = uniform; Ju et al. 1995)
        crosstalk: spectral crosstalk fraction between channels
                   (0 = none, 0.03 = realistic; Huang et al. 1992)
        baseline_drift: sinusoidal baseline drift amplitude
                        (0 = none, 15 = mild wandering)
        smooth: post-generation Gaussian smoothing sigma
                (0 = none, 1.0 = mild; Virtanen et al. 2020)

    Returns:
        traces: dict {"A","C","G","T" → numpy int16 arrays}
        base_calls: str of called bases
        peak_positions: list of int peak center positions
        quality_scores: list of int quality values (phred scale)
    """
    n = len(base_freq)
    # Cap spacing so peak positions fit in uint16 (AB1 PLOC format limit)
    max_uint16 = 65535
    if n * spacing > max_uint16:
        spacing = max(1, max_uint16 // n)
    signal_len = n * spacing

    # Wider window for asymmetric peaks
    half_window = int(3 * sigma * max(asymmetry, 1.0)) + 1

    # Initialize float traces
    trace = {b: np.zeros(signal_len, dtype=np.float64) for b in BASES}

    # --- Phase 1.1: Generate peaks (symmetric or asymmetric) ---
    use_asymmetric = (asymmetry != 1.0)
    for i, freq in enumerate(base_freq):
        center = i * spacing
        lo = max(0, center - half_window)
        hi = min(signal_len, center + half_window + 1)
        x = np.arange(lo, hi)

        if use_asymmetric:
            g = _asymmetric_gaussian(x, center, sigma, asymmetry)
        else:
            g = _gaussian(x, center, sigma)

        for base in BASES:
            height = freq.get(base, 0.0) * scale
            trace[base][lo:hi] += height * g

    # --- Phase 1.2: Dye-specific amplitude scaling ---
    if dye_scaling is not None:
        for base in BASES:
            factor = dye_scaling.get(base, 1.0)
            trace[base] *= factor

    # --- Phase 1.3: Cross-channel spectral crosstalk ---
    if crosstalk > 0:
        # Each channel bleeds a fraction into every other channel
        original = {b: trace[b].copy() for b in BASES}
        for b1 in BASES:
            for b2 in BASES:
                if b1 != b2:
                    trace[b1] += crosstalk * original[b2]

    # --- Phasing (tailing effect) ---
    if phasing > 0:
        for base in BASES:
            for i in range(1, signal_len):
                trace[base][i] += phasing * trace[base][i - 1]

    # --- Noise ---
    if noise > 0:
        rng = np.random.default_rng(42)
        for base in BASES:
            trace[base] += rng.normal(0, noise, signal_len)

    # --- Signal decay ---
    if decay > 0:
        decay_curve = np.exp(-np.linspace(0, decay, signal_len))
        for base in BASES:
            trace[base] *= decay_curve

    # --- Phase 1.4: Baseline drift ---
    if baseline_drift > 0:
        rng_drift = np.random.default_rng(123)
        x_norm = np.linspace(0, 1, signal_len)
        for base in BASES:
            # Low-frequency sinusoidal drift (2-5 cycles) with random phase
            n_cycles = rng_drift.uniform(2, 5)
            phase = rng_drift.uniform(0, 2 * np.pi)
            drift = baseline_drift * np.sin(2 * np.pi * n_cycles * x_norm + phase)
            trace[base] += drift

    # --- Phase 1.5: Local smoothing ---
    if smooth > 0:
        try:
            from scipy.ndimage import gaussian_filter1d
            for base in BASES:
                trace[base] = gaussian_filter1d(trace[base], sigma=smooth)
        except ImportError:
            pass  # scipy not available, skip smoothing

    # --- Clamp to int16 range ---
    for base in BASES:
        trace[base] = np.clip(trace[base], 0, 32767).astype(np.int16)

    # --- Base calls, peak positions, and quality scores ---
    base_calls = ""
    peak_positions = []
    quality_scores = []
    for i, freq in enumerate(base_freq):
        called = max(BASES, key=lambda b: freq.get(b, 0.0))
        base_calls += called
        peak_positions.append(i * spacing)

        # Quality = -10 * log10(1 - freq_of_called_base), capped at 60
        called_freq = freq.get(called, 0.0)
        if called_freq >= 1.0:
            q = 60
        elif called_freq <= 0.0:
            q = 0
        else:
            q = min(60, int(-10 * math.log10(1 - called_freq)))
        quality_scores.append(q)

    return trace, base_calls, peak_positions, quality_scores
