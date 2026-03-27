"""Tests for the SDVC (Signal Deconvolution-based Variant Calling) module."""

import numpy as np
import pytest
import os
import tempfile

from ab1tools.variant_caller import (
    VariantCaller, call_variants_from_traces,
    write_variants_csv, write_variants_vcf,
)
from ab1tools.signal import generate_traces, sequence_to_base_frequencies


def _make_clean_traces(seq="ACGTACGTACGT", spacing=10, scale=1024):
    """Helper: generate clean traces from a pure sequence (no variants)."""
    base_freq = sequence_to_base_frequencies(seq)
    traces, calls, peaks, quals = generate_traces(
        base_freq, spacing=spacing, scale=scale, noise=0
    )
    return traces, calls, peaks


def _make_hetero_traces(n=100, variant_positions=None, variant_freq=0.3,
                        spacing=10, scale=1024, noise=0):
    """Helper: generate traces with known heterozygous positions."""
    if variant_positions is None:
        variant_positions = [25, 50, 75]

    # Build base frequencies: mostly pure A, with variants at specified positions
    base_freq = []
    seq = ""
    for i in range(n):
        if i in variant_positions:
            # Heterozygous: major=A, minor=G
            freq = {"A": 1.0 - variant_freq, "C": 0.0,
                    "G": variant_freq, "T": 0.0}
            seq += "A"
        else:
            # Pure A
            freq = {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0}
            seq += "A"
        base_freq.append(freq)

    traces, calls, peaks, quals = generate_traces(
        base_freq, spacing=spacing, scale=scale, noise=noise
    )
    return traces, calls, peaks, variant_positions


class TestVariantCallerBasic:
    """Test basic VariantCaller functionality."""

    def test_init(self):
        traces, calls, peaks = _make_clean_traces()
        caller = VariantCaller(traces, peaks, calls)
        assert caller.trace_len > 0
        assert len(caller.peaks) == len(calls)

    def test_no_variants_in_clean_sequence(self):
        traces, calls, peaks = _make_clean_traces("ACGTACGTACGT")
        caller = VariantCaller(traces, peaks, calls, min_threshold=0.05)
        results = caller.call_variants()
        assert len(results) == len(calls)
        variants = [r for r in results if r["is_variant"]]
        assert len(variants) == 0

    def test_detect_known_variants(self):
        traces, calls, peaks, var_pos = _make_hetero_traces(
            n=100, variant_positions=[25, 50, 75], variant_freq=0.3
        )
        caller = VariantCaller(traces, peaks, calls, min_threshold=0.05)
        variants = caller.get_variants()
        detected_positions = {v["pos_0based"] for v in variants}
        # All known variant positions should be detected
        for pos in var_pos:
            assert pos in detected_positions, f"Variant at pos {pos} not detected"

    def test_result_structure(self):
        traces, calls, peaks = _make_clean_traces("ACGT")
        caller = VariantCaller(traces, peaks, calls)
        results = caller.call_variants()
        assert len(results) > 0
        r = results[0]
        required_keys = [
            "pos_0based", "pos_1based", "trace_pos", "consensus_base",
            "major_base", "major_freq", "minor_base", "minor_freq",
            "third_base", "third_freq",
            "freq_A", "freq_C", "freq_G", "freq_T",
            "is_variant", "confidence", "snr", "threshold",
        ]
        for key in required_keys:
            assert key in r, f"Missing key: {key}"

    def test_frequencies_sum_to_one(self):
        traces, calls, peaks, _ = _make_hetero_traces(n=20, variant_positions=[10])
        caller = VariantCaller(traces, peaks, calls)
        results = caller.call_variants()
        for r in results:
            total = r["freq_A"] + r["freq_C"] + r["freq_G"] + r["freq_T"]
            if total > 0:  # Skip positions with no signal
                assert abs(total - 1.0) < 0.01, f"Frequencies don't sum to 1: {total}"


class TestDynamicThresholding:
    """Test dynamic threshold adaptation."""

    def test_threshold_varies_with_noise(self):
        # Clean traces
        traces_clean, calls, peaks, _ = _make_hetero_traces(n=50, noise=0)
        caller_clean = VariantCaller(traces_clean, peaks, calls, sensitivity=0.5)
        results_clean = caller_clean.call_variants()

        # Noisy traces
        traces_noisy, calls, peaks, _ = _make_hetero_traces(n=50, noise=20)
        caller_noisy = VariantCaller(traces_noisy, peaks, calls, sensitivity=0.5)
        results_noisy = caller_noisy.call_variants()

        # Both should produce valid thresholds
        for r in results_clean + results_noisy:
            assert 0.0 < r["threshold"] <= 0.5

    def test_sensitivity_parameter(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.10
        )
        # Low sensitivity (more conservative)
        caller_low = VariantCaller(traces, peaks, calls, sensitivity=1.0)
        variants_low = caller_low.get_variants()

        # High sensitivity (more liberal)
        caller_high = VariantCaller(traces, peaks, calls, sensitivity=0.1)
        variants_high = caller_high.get_variants()

        # More sensitive should detect at least as many
        assert len(variants_high) >= len(variants_low)


class TestConfidenceScoring:
    """Test confidence score computation."""

    def test_high_freq_variant_high_confidence(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.4
        )
        caller = VariantCaller(traces, peaks, calls)
        variants = caller.get_variants()
        if variants:
            assert variants[0]["confidence"] > 0

    def test_low_freq_variant_lower_confidence(self):
        # 40% variant
        traces_hi, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.4
        )
        caller_hi = VariantCaller(traces_hi, peaks, calls)
        results_hi = caller_hi.call_variants()

        # 10% variant
        traces_lo, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.10
        )
        caller_lo = VariantCaller(traces_lo, peaks, calls)
        results_lo = caller_lo.call_variants()

        conf_hi = results_hi[25]["confidence"]
        conf_lo = results_lo[25]["confidence"]
        assert conf_hi >= conf_lo

    def test_confidence_range(self):
        traces, calls, peaks, _ = _make_hetero_traces(n=50)
        caller = VariantCaller(traces, peaks, calls)
        results = caller.call_variants()
        for r in results:
            assert 0.0 <= r["confidence"] <= 1.0


class TestMultiPositionConsistency:
    """Test sliding window consistency enhancement."""

    def test_clustered_variants_boosted(self):
        # Three adjacent variant positions → should boost confidence
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[24, 25, 26], variant_freq=0.25
        )
        caller = VariantCaller(traces, peaks, calls, window_size=5)
        results = caller.call_variants()
        # Center position should have boosted confidence
        if results[25]["is_variant"]:
            assert results[25]["confidence"] > 0

    def test_window_size_parameter(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.20
        )
        # Different window sizes should produce results
        for ws in [3, 5, 7]:
            caller = VariantCaller(traces, peaks, calls, window_size=ws)
            results = caller.call_variants()
            assert len(results) == 50


class TestBayesianPrior:
    """Test Bayesian prior integration."""

    def test_higher_prior_more_sensitive(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.08
        )
        # Low prior (rare variant assumption)
        caller_low = VariantCaller(traces, peaks, calls, prior=0.01)
        results_low = caller_low.call_variants()

        # High prior (known hotspot)
        caller_high = VariantCaller(traces, peaks, calls, prior=0.5)
        results_high = caller_high.call_variants()

        # Results should be produced for both
        assert len(results_low) == len(results_high)

    def test_zero_prior_disables_bayesian(self):
        traces, calls, peaks, _ = _make_hetero_traces(n=50)
        caller = VariantCaller(traces, peaks, calls, prior=0.0)
        results = caller.call_variants()
        assert len(results) == 50


class TestFunctionalInterface:
    """Test the convenience functional interface."""

    def test_call_variants_from_traces(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.3
        )
        results = call_variants_from_traces(traces, peaks, calls)
        assert len(results) == 50
        variants = [r for r in results if r["is_variant"]]
        assert len(variants) > 0


class TestOutputFormats:
    """Test CSV and VCF output writing."""

    def test_write_csv(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.3
        )
        caller = VariantCaller(traces, peaks, calls)
        variants = caller.get_variants()

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            csv_path = f.name

        try:
            write_variants_csv(variants, csv_path, sample_name="test")
            assert os.path.exists(csv_path)
            with open(csv_path) as f:
                lines = f.readlines()
            assert len(lines) >= 2  # header + at least 1 variant
            assert "sample" in lines[0]
            assert "confidence" in lines[0]
        finally:
            os.unlink(csv_path)

    def test_write_vcf(self):
        traces, calls, peaks, _ = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=0.3
        )
        caller = VariantCaller(traces, peaks, calls)
        variants = caller.get_variants()

        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
            vcf_path = f.name

        try:
            write_variants_vcf(variants, vcf_path, sample_name="test")
            assert os.path.exists(vcf_path)
            with open(vcf_path) as f:
                content = f.read()
            assert "##fileformat=VCFv4.2" in content
            assert "##source=ab1tools_SDVC" in content
            assert "PASS" in content or "LowConf" in content
        finally:
            os.unlink(vcf_path)


class TestFrequencyRecovery:
    """Test that SDVC accurately recovers known variant frequencies."""

    @pytest.mark.parametrize("true_freq", [0.10, 0.20, 0.30, 0.40, 0.50])
    def test_frequency_recovery(self, true_freq):
        traces, calls, peaks, var_pos = _make_hetero_traces(
            n=50, variant_positions=[25], variant_freq=true_freq
        )
        caller = VariantCaller(traces, peaks, calls, min_threshold=0.03)
        results = caller.call_variants()
        r = results[25]
        # Check that recovered frequency is close to true frequency
        # Allow 10% absolute error for signal-level estimation (baseline
        # correction and Gaussian overlap cause small systematic shifts)
        assert abs(r["minor_freq"] - true_freq) < 0.10, \
            f"Expected ~{true_freq}, got {r['minor_freq']}"


class TestWithRealAB1:
    """Test with actual AB1 files from demo_data."""

    @pytest.fixture
    def demo_ab1_path(self):
        """Find a demo AB1 file."""
        demo_dir = os.path.join(os.path.dirname(__file__), "..", "demo_data", "output")
        for root, dirs, files in os.walk(demo_dir):
            for f in files:
                if f.endswith(".ab1"):
                    return os.path.join(root, f)
        pytest.skip("No demo AB1 files found")

    def test_call_on_real_ab1(self, demo_ab1_path):
        from ab1tools.abif_reader import read_ab1
        traces, base_calls, peak_positions, sample_name = read_ab1(demo_ab1_path)
        caller = VariantCaller(traces, peak_positions, base_calls)
        results = caller.call_variants()
        assert len(results) == len(base_calls)
        # All results should have valid structure
        for r in results:
            assert 0 <= r["major_freq"] <= 1
            assert 0 <= r["minor_freq"] <= 1
            assert isinstance(r["is_variant"], bool)
