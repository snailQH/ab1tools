"""ABIF (AB1) binary file reader.

Reads trace data, base calls, and peak positions from AB1 files
for re-plotting as PNG chromatograms.
"""

import struct
import numpy as np

ABIF_MAGIC = b"ABIF"
DIR_ENTRY_SIZE = 28


def _read_dir_entries(f, dir_offset, num_entries):
    """Read directory entries from an AB1 file."""
    f.seek(dir_offset)
    entries = {}
    for _ in range(num_entries):
        raw = f.read(DIR_ENTRY_SIZE)
        tag_name = raw[0:4].decode("ascii")
        tag_num = struct.unpack(">I", raw[4:8])[0]
        elem_type = struct.unpack(">H", raw[8:10])[0]
        elem_size = struct.unpack(">H", raw[10:12])[0]
        num_elem = struct.unpack(">I", raw[12:16])[0]
        data_size = struct.unpack(">I", raw[16:20])[0]
        data_offset_raw = raw[20:24]
        # dataHandle = raw[24:28]  # unused

        key = f"{tag_name}{tag_num}"

        if data_size <= 4:
            # Data is inline in the offset field
            data = data_offset_raw[:data_size]
        else:
            data_offset = struct.unpack(">I", data_offset_raw)[0]
            pos = f.tell()
            f.seek(data_offset)
            data = f.read(data_size)
            f.seek(pos)

        entries[key] = {
            "tag": tag_name,
            "num": tag_num,
            "type": elem_type,
            "elem_size": elem_size,
            "num_elem": num_elem,
            "data": data,
        }

    return entries


def read_ab1(filename):
    """Read an AB1 file and extract traces, base calls, and peak positions.

    Returns:
        traces: dict with keys "A","C","G","T" -> numpy int16 arrays
        base_calls: str of base calls
        peak_positions: list of int peak positions
        sample_name: str sample name
    """
    with open(filename, "rb") as f:
        # Read header
        magic = f.read(4)
        if magic != ABIF_MAGIC:
            raise ValueError(f"Not an ABIF file: {filename}")

        version = struct.unpack(">H", f.read(2))[0]

        # Read root directory entry
        root_raw = f.read(DIR_ENTRY_SIZE)
        # root tag name = "tdir", number = 1
        num_entries = struct.unpack(">I", root_raw[12:16])[0]
        dir_offset = struct.unpack(">I", root_raw[20:24])[0]

        entries = _read_dir_entries(f, dir_offset, num_entries)

    # Extract trace channels (DATA9=G, DATA10=A, DATA11=T, DATA12=C)
    traces = {}
    channel_map = {"DATA9": "G", "DATA10": "A", "DATA11": "T", "DATA12": "C"}
    for tag_key, base in channel_map.items():
        if tag_key in entries:
            data = entries[tag_key]["data"]
            traces[base] = np.frombuffer(data, dtype=">i2").astype(np.int16)

    if not traces:
        raise ValueError(f"No trace data (DATA9-12) found in {filename}")

    # Extract base calls (try PBAS1, then PBAS2)
    base_calls = ""
    for key in ["PBAS1", "PBAS2"]:
        if key in entries:
            base_calls = entries[key]["data"].decode("ascii")
            break

    # Extract peak positions (try PLOC1, then PLOC2)
    peak_positions = []
    for key in ["PLOC1", "PLOC2"]:
        if key in entries:
            data = entries[key]["data"]
            peak_positions = np.frombuffer(data, dtype=">u2").astype(int).tolist()
            break

    # Extract sample name
    sample_name = "Sample"
    if "SMPL1" in entries:
        raw = entries["SMPL1"]["data"]
        # Pascal string: first byte is length
        if raw:
            str_len = raw[0]
            sample_name = raw[1:1 + str_len].decode("ascii", errors="replace")

    return traces, base_calls, peak_positions, sample_name
