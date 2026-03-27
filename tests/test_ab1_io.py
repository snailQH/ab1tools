"""Tests for AB1 read/write round-trip and format conversion."""
import os
import tempfile
import numpy as np
import pytest
from ab1tools.signal import generate_traces, sequence_to_base_frequencies
from ab1tools.abif_writer import write_ab1
from ab1tools.abif_reader import read_ab1


class TestAB1RoundTrip:
    def test_write_read_roundtrip(self):
        bf = sequence_to_base_frequencies("ACGTACGTACGT")
        traces, calls, peaks, quals = generate_traces(bf)
        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            path = f.name
        try:
            write_ab1(path, traces, calls, peaks, sample_name="test", quality_scores=quals)
            t2, c2, p2, name = read_ab1(path)
            assert c2 == calls
            assert name == "test"
            assert len(p2) == len(peaks)
            for b in "ACGT":
                np.testing.assert_array_equal(t2[b], traces[b])
        finally:
            os.unlink(path)

    def test_abif_magic(self):
        bf = sequence_to_base_frequencies("ACGT")
        traces, calls, peaks, quals = generate_traces(bf)
        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            path = f.name
        try:
            write_ab1(path, traces, calls, peaks)
            with open(path, "rb") as f:
                magic = f.read(4)
            assert magic == b"ABIF"
        finally:
            os.unlink(path)

    def test_long_sequence(self):
        bf = sequence_to_base_frequencies("A" * 5000)
        traces, calls, peaks, quals = generate_traces(bf)
        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            path = f.name
        try:
            write_ab1(path, traces, calls, peaks, quality_scores=quals)
            t2, c2, p2, _ = read_ab1(path)
            assert len(c2) == 5000
        finally:
            os.unlink(path)

    def test_realistic_roundtrip(self):
        bf = sequence_to_base_frequencies("ACGTACGT")
        traces, calls, peaks, quals = generate_traces(
            bf, noise=5, phasing=0.1, decay=0.5, asymmetry=1.3,
            crosstalk=0.03, baseline_drift=10, smooth=1.0,
        )
        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            path = f.name
        try:
            write_ab1(path, traces, calls, peaks, sample_name="realistic",
                      quality_scores=quals)
            t2, c2, p2, name = read_ab1(path)
            assert c2 == calls
            assert name == "realistic"
        finally:
            os.unlink(path)

    def test_invalid_file_raises(self):
        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            f.write(b"NOT_ABIF_DATA")
            path = f.name
        try:
            with pytest.raises(ValueError, match="Not an ABIF"):
                read_ab1(path)
        finally:
            os.unlink(path)


class TestDemoAB1Files:
    """Test with actual demo data AB1 files."""

    @pytest.fixture
    def demo_ab1_files(self):
        demo_dir = os.path.join(os.path.dirname(__file__), "..", "demo_data", "output")
        files = []
        for root, dirs, fnames in os.walk(demo_dir):
            for f in fnames:
                if f.endswith(".ab1"):
                    files.append(os.path.join(root, f))
        if not files:
            pytest.skip("No demo AB1 files")
        return files[:5]  # test first 5

    def test_all_demo_files_readable(self, demo_ab1_files):
        for path in demo_ab1_files:
            traces, calls, peaks, name = read_ab1(path)
            assert len(calls) > 0
            assert len(peaks) == len(calls)
            for b in "ACGT":
                assert b in traces
                assert len(traces[b]) > 0
