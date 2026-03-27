"""Tests for signal generation module (Phase 1 realism features)."""
import numpy as np
import pytest
from ab1tools.signal import (
    generate_traces, sequence_to_base_frequencies,
    _gaussian, _asymmetric_gaussian, compute_snr,
)


class TestGaussianFunctions:
    def test_symmetric_gaussian_peak(self):
        x = np.arange(0, 20)
        g = _gaussian(x, 10, 2.0)
        assert g[10] == pytest.approx(1.0)
        assert g[0] < 0.01

    def test_asymmetric_gaussian_symmetric_when_1(self):
        x = np.arange(0, 20)
        sym = _gaussian(x, 10, 2.0)
        asym = _asymmetric_gaussian(x, 10, 2.0, 1.0)
        np.testing.assert_array_almost_equal(sym, asym)

    def test_asymmetric_gaussian_rightward_tail(self):
        x = np.arange(0, 20)
        asym = _asymmetric_gaussian(x, 10, 2.0, 1.3)
        # Right side should be broader (higher values further from center)
        assert asym[13] > _gaussian(np.array([13]), 10, 2.0)[0]
        # Left side unchanged
        assert asym[7] == pytest.approx(_gaussian(np.array([7]), 10, 2.0)[0])


class TestGenerateTraces:
    def test_basic_output(self):
        bf = sequence_to_base_frequencies("ACGT")
        traces, calls, peaks, quals = generate_traces(bf)
        assert len(calls) == 4
        assert calls == "ACGT"
        assert len(peaks) == 4
        assert "A" in traces and "C" in traces

    def test_spacing_overflow_protection(self):
        bf = sequence_to_base_frequencies("A" * 7000)
        traces, calls, peaks, quals = generate_traces(bf, spacing=10)
        assert max(peaks) <= 65535

    def test_asymmetry_parameter(self):
        bf = sequence_to_base_frequencies("AAAA")
        t_sym, _, _, _ = generate_traces(bf, asymmetry=1.0)
        t_asym, _, _, _ = generate_traces(bf, asymmetry=1.3)
        # Asymmetric should differ from symmetric
        assert not np.array_equal(t_sym["A"], t_asym["A"])

    def test_dye_scaling(self):
        bf = [{"A":1,"C":0,"G":0,"T":0}, {"A":0,"C":1,"G":0,"T":0},
              {"A":0,"C":0,"G":1,"T":0}, {"A":0,"C":0,"G":0,"T":1}]
        t_unscaled, _, _, _ = generate_traces(bf, dye_scaling=None)
        t_scaled, _, _, _ = generate_traces(bf, dye_scaling={"G":0.5,"A":1.0,"T":1.0,"C":1.0})
        # G channel should be halved
        assert np.max(t_scaled["G"]) < np.max(t_unscaled["G"])

    def test_crosstalk(self):
        bf = [{"A":1,"C":0,"G":0,"T":0}] * 5
        t_no, _, _, _ = generate_traces(bf, crosstalk=0.0)
        t_ct, _, _, _ = generate_traces(bf, crosstalk=0.1)
        # With crosstalk, non-A channels should have signal
        assert np.max(t_no["C"]) == 0
        assert np.max(t_ct["C"]) > 0

    def test_baseline_drift(self):
        bf = sequence_to_base_frequencies("A" * 20)
        t_no, _, _, _ = generate_traces(bf, baseline_drift=0)
        t_drift, _, _, _ = generate_traces(bf, baseline_drift=50)
        assert not np.array_equal(t_no["A"], t_drift["A"])

    def test_smoothing(self):
        bf = sequence_to_base_frequencies("ACGT" * 5)
        t_no, _, _, _ = generate_traces(bf, smooth=0)
        t_sm, _, _, _ = generate_traces(bf, smooth=2.0)
        # Smoothed should have lower max (averaging reduces peaks)
        assert np.max(t_sm["A"]) <= np.max(t_no["A"])

    def test_noise(self):
        bf = sequence_to_base_frequencies("A" * 10)
        t1, _, _, _ = generate_traces(bf, noise=0)
        t2, _, _, _ = generate_traces(bf, noise=10)
        # With noise, inter-peak regions should have non-zero values
        assert np.std(t2["C"].astype(float)) > np.std(t1["C"].astype(float))

    def test_decay(self):
        bf = sequence_to_base_frequencies("A" * 50)
        traces, _, peaks, _ = generate_traces(bf, decay=2.0)
        # First peak should be taller than last
        first = traces["A"][peaks[0]]
        last = traces["A"][peaks[-1]]
        assert first > last

    def test_phasing(self):
        bf = sequence_to_base_frequencies("ACGT")
        t_no, _, _, _ = generate_traces(bf, phasing=0)
        t_ph, _, _, _ = generate_traces(bf, phasing=0.2)
        assert not np.array_equal(t_no["A"], t_ph["A"])

    def test_full_realistic_preset(self):
        bf = sequence_to_base_frequencies("ACGTACGT")
        traces, calls, peaks, quals = generate_traces(
            bf, noise=5, phasing=0.1, decay=0.5,
            asymmetry=1.3, dye_scaling={"G":0.8,"A":1.0,"T":0.95,"C":0.85},
            crosstalk=0.03, baseline_drift=15, smooth=1.0,
        )
        assert calls == "ACGTACGT"
        assert len(peaks) == 8
        for b in "ACGT":
            assert traces[b].dtype == np.int16


class TestComputeSNR:
    def test_clean_trace_high_snr(self):
        bf = sequence_to_base_frequencies("ACGT" * 10)
        traces, _, _, _ = generate_traces(bf, noise=0)
        snr = compute_snr(traces)
        for b in "ACGT":
            assert snr[b] > 100  # High SNR for clean traces

    def test_noisy_trace_lower_snr(self):
        bf = sequence_to_base_frequencies("ACGT" * 10)
        t_clean, _, _, _ = generate_traces(bf, noise=0)
        t_noisy, _, _, _ = generate_traces(bf, noise=20)
        snr_clean = compute_snr(t_clean)
        snr_noisy = compute_snr(t_noisy)
        # At least some channels should have lower SNR with noise
        assert min(snr_noisy.values()) <= max(snr_clean.values())
