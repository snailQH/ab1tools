"""Chromatogram PNG visualization."""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Standard Sanger trace colors
COLORS = {"A": "#00CC00", "C": "#0000FF", "G": "#000000", "T": "#FF0000"}

# Bases per row in the detail PNG
BASES_PER_ROW = 100


def _base_position_formatter(spacing):
    """Return a matplotlib FuncFormatter that converts trace position to base position."""
    def formatter(x, pos):
        base_pos = int(round(x / spacing)) + 1  # 1-based
        return str(base_pos)
    return FuncFormatter(formatter)


def _plot_trace_panel(ax, traces, base_calls, peak_positions, spacing,
                      xlo, xhi, title, show_bases=True, linewidth=0.8,
                      base_fontsize=5, show_pos_numbers=False):
    """Plot a single trace panel on the given axes."""
    x = np.arange(xlo, xhi)

    for base in "ACGT":
        ax.plot(x, traces[base][xlo:xhi], color=COLORS[base], label=base,
                linewidth=linewidth, alpha=0.85)

    if show_bases:
        ymax = max(traces[b][xlo:xhi].max() for b in "ACGT")
        if ymax <= 0:
            ymax = 1
        y_text = ymax * 1.08
        y_pos_num = -ymax * 0.08

        for i, pos in enumerate(peak_positions):
            if xlo <= pos < xhi and i < len(base_calls):
                base = base_calls[i]
                ax.text(pos, y_text, base, ha="center", va="bottom",
                        fontsize=base_fontsize, color=COLORS.get(base, "gray"),
                        fontweight="bold", fontfamily="monospace")

                # Show base position number every 10 bases
                if show_pos_numbers and (i + 1) % 10 == 0:
                    ax.text(pos, y_pos_num, str(i + 1), ha="center", va="top",
                            fontsize=4, color="gray", fontfamily="monospace")

    ax.set_xlim(xlo, xhi)
    ax.set_ylabel("Signal", fontsize=8)
    ax.set_title(title, fontsize=9, pad=15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Show base positions on x-axis instead of trace positions
    ax.xaxis.set_major_formatter(_base_position_formatter(spacing))


def plot_chromatogram(traces, base_calls, peak_positions, output_path,
                      title="Chromatogram", window=None, figsize=(20, 5),
                      show_bases=True, dpi=150, bases_per_row=BASES_PER_ROW,
                      base_range=None):
    """Generate chromatogram PNGs covering a sequence region.

    Args:
        base_range: optional (start_0based, end_0based) inclusive tuple to
                    restrict output to a specific region. None = full sequence.

    Produces a multi-page PNG:
      - Each page has up to 4 rows of bases_per_row bases each
      - Base letters and position numbers annotated on every row
      - Full sequence (or specified region) is covered across PNG files

    Output files: sample.png, sample_p2.png, sample_p3.png, etc.
    """
    trace_len = len(traces["A"])
    n_bases = len(base_calls)
    spacing = peak_positions[1] - peak_positions[0] if len(peak_positions) > 1 else 10

    # Apply base_range filter
    if base_range is not None:
        range_start, range_end = base_range
        range_start = max(0, range_start)
        range_end = min(n_bases - 1, range_end)
    else:
        range_start = 0
        range_end = n_bases - 1

    region_bases = range_end - range_start + 1

    rows_per_page = 4
    total_rows = (region_bases + bases_per_row - 1) // bases_per_row
    total_pages = (total_rows + rows_per_page - 1) // rows_per_page

    output_paths = []
    base_name, ext = os.path.splitext(output_path)

    for page in range(total_pages):
        start_row = page * rows_per_page
        end_row = min(start_row + rows_per_page, total_rows)
        n_rows = end_row - start_row

        fig, axes = plt.subplots(n_rows, 1, figsize=(20, 3 * n_rows),
                                 squeeze=False)
        page_label = f" (page {page + 1}/{total_pages})" if total_pages > 1 else ""
        fig.suptitle(f"{title}{page_label}", fontsize=13, fontweight="bold", y=0.99)

        for row_idx in range(n_rows):
            ax = axes[row_idx, 0]
            global_row = start_row + row_idx
            base_start = range_start + global_row * bases_per_row
            base_end = min(base_start + bases_per_row, range_end + 1)

            xlo = base_start * spacing
            xhi = min(base_end * spacing, trace_len)

            row_title = f"Bases {base_start + 1}–{base_end}"
            _plot_trace_panel(ax, traces, base_calls, peak_positions, spacing,
                              xlo, xhi, row_title,
                              show_bases=True, linewidth=1.0,
                              base_fontsize=6, show_pos_numbers=True)

            if row_idx == 0:
                ax.legend(loc="upper right", fontsize=7)

        if n_rows > 0:
            axes[-1, 0].set_xlabel("Base Position", fontsize=9)

        plt.tight_layout(rect=[0, 0, 1, 0.97])

        if page == 0:
            page_path = output_path
        else:
            page_path = f"{base_name}_p{page + 1}{ext}"

        fig.savefig(page_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        output_paths.append(page_path)

    return output_paths


def plot_hetero_windows(traces, base_calls, peak_positions, hetero_sites,
                        windows, output_path, title="Chromatogram",
                        dpi=150):
    """Generate a single PNG showing only the regions around heterozygous sites.

    Each window gets one row. Hetero positions are highlighted with a
    vertical shaded bar.
    """
    if not windows:
        return None

    trace_len = len(traces["A"])
    spacing = peak_positions[1] - peak_positions[0] if len(peak_positions) > 1 else 10
    n_rows = len(windows)

    fig, axes = plt.subplots(n_rows, 1, figsize=(20, 3.5 * n_rows), squeeze=False)
    fig.suptitle(f"{title} — Heterozygous Sites", fontsize=13, fontweight="bold", y=0.99)

    for row_idx, (base_start, base_end) in enumerate(windows):
        ax = axes[row_idx, 0]
        xlo = base_start * spacing
        xhi = min((base_end + 1) * spacing, trace_len)

        # Find hetero sites in this window for the title
        sites_in_window = [s for s in hetero_sites
                           if base_start <= s["pos_0based"] <= base_end]
        pos_strs = []
        for s in sites_in_window:
            pos_strs.append(
                f"pos {s['pos_1based']}: {s['major_base']}={s['major_freq']*100:.0f}%"
                f"/{s['minor_base']}={s['minor_freq']*100:.0f}%"
            )
        row_title = f"Bases {base_start + 1}–{base_end + 1}  |  {', '.join(pos_strs)}"

        _plot_trace_panel(ax, traces, base_calls, peak_positions, spacing,
                          xlo, xhi, row_title,
                          show_bases=True, linewidth=1.2,
                          base_fontsize=7, show_pos_numbers=True)

        # Highlight hetero positions with vertical shading
        for s in sites_in_window:
            center = s["pos_0based"] * spacing
            ax.axvspan(center - spacing * 0.5, center + spacing * 0.5,
                       alpha=0.15, color="orange", zorder=0)

        if row_idx == 0:
            ax.legend(loc="upper right", fontsize=7)

    if n_rows > 0:
        axes[-1, 0].set_xlabel("Base Position", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path
