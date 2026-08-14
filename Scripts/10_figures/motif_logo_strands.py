#!/usr/bin/env python3

"""
Strand-aware DNA motif visualization using Logomaker.

The input DNA sequence is assumed to be written in the positive-strand
orientation.

The figure displays:

- The complete positive-strand sequence as small letters above zero.
- The complete complementary negative strand as small inverted letters
  below zero.
- Positive-strand motif hits enlarged above zero.
- Negative-strand motif hits enlarged below zero.
- All motif matches above a user-defined PWM score threshold.
- SVG, PDF, and PNG output.

Run from macOS Terminal:

    python3 motif_logo_strands.py
"""

# ============================================================
# INSTALL MISSING PACKAGES
# ============================================================

import sys
import subprocess
import importlib.util


def install_if_missing(import_name, pip_name=None):
    """
    Install a Python package if it is not already available.
    """

    if pip_name is None:
        pip_name = import_name

    if importlib.util.find_spec(import_name) is None:
        print(f"Installing missing package: {pip_name}")

        subprocess.check_call(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--user",
                pip_name
            ]
        )


for package_name in [
    "numpy",
    "pandas",
    "matplotlib",
    "logomaker"
]:
    install_if_missing(package_name)

# ============================================================
# IMPORTS
# ============================================================

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

# ============================================================
# USER INPUT
# ============================================================

# Input sequence written in the positive-strand orientation.
sequence = (
    "TTTTCCTGACACAATTCAGT"
)

# Background nucleotide frequencies.
background = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25
}

# ============================================================
# MOTIF THRESHOLD
# ============================================================

# Hits are retained when:
#
# hit score >= score_threshold_fraction × theoretical maximum score
#
# Suggested values:
#
# 0.40 = permissive
# 0.50 = moderate
# 0.60 = moderately stringent
# 0.70 = stringent
# 0.80 = very stringent

score_threshold_fraction = 0.70

# Optional maximum number of hits.
#
# None:
#     Retain all hits above the threshold.
#
# Integer:
#     Retain only the highest-scoring N hits.
#
# Example:
#     top_n_hits = 5

top_n_hits = None

# Scan the reverse-complement strand.
scan_reverse_strand = True

# ============================================================
# PWM
# ============================================================

# Rows correspond to motif positions.
# Columns correspond to A, C, G, T probabilities.

pwm = pd.DataFrame(
    [
        [0.041, 0.029, 0.029, 0.900],
        [0.037, 0.026, 0.900, 0.037],
        [0.900, 0.029, 0.029, 0.041],
        [0.037, 0.900, 0.026, 0.037],
        [0.291, 0.209, 0.209, 0.291],
        [0.291, 0.209, 0.209, 0.291],
        [0.900, 0.029, 0.029, 0.041],
        [0.058, 0.450, 0.042, 0.450],
    ],
    columns=list("ACGT")
)

# ============================================================
# PLOT SETTINGS
# ============================================================

# Height of bases that are not part of a retained motif match.
baseline_height = 0.16

# Minimum and maximum height of retained motif bases.
motif_min_height = 0.35
motif_max_height = 2.50

# Transparency of highlighted motif regions.
highlight_alpha = 0.11

# Highlight colors.
plus_highlight_color = "gold"
minus_highlight_color = "lightskyblue"

# Figure dimensions.
#
# Width is also adjusted automatically according to sequence length.
minimum_figure_width = 10
width_per_base = 0.32
figure_height = 6.0

# Export filenames.
svg_output = "motif_logo_strands.svg"
pdf_output = "motif_logo_strands.pdf"
png_output = "motif_logo_strands.png"

# Show the plot in a macOS window after saving.
show_plot = True

# ============================================================
# VALIDATION
# ============================================================


def clean_and_validate_inputs(sequence, pwm):
    """
    Clean and validate the DNA sequence and PWM.
    """

    clean_sequence = (
        sequence.upper()
        .replace(" ", "")
        .replace("\n", "")
        .replace("\t", "")
        .replace("\r", "")
    )

    if not clean_sequence:
        raise ValueError("The DNA sequence is empty.")

    invalid_bases = sorted(
        set(clean_sequence) - set("ACGT")
    )

    if invalid_bases:
        raise ValueError(
            "The DNA sequence contains unsupported characters: "
            + ", ".join(invalid_bases)
        )

    if len(clean_sequence) < len(pwm):
        raise ValueError(
            f"The DNA sequence is {len(clean_sequence)} bp long, "
            f"but the PWM is {len(pwm)} bp long."
        )

    if list(pwm.columns) != list("ACGT"):
        raise ValueError(
            "PWM columns must be ordered as A, C, G, T."
        )

    if (pwm < 0).any().any():
        raise ValueError(
            "PWM probabilities cannot be negative."
        )

    pwm_row_sums = pwm.sum(axis=1)

    if not np.allclose(
        pwm_row_sums,
        1.0,
        atol=0.01
    ):
        raise ValueError(
            "Each PWM row must sum to approximately 1. "
            f"Current row sums: {pwm_row_sums.tolist()}"
        )

    if not 0 <= score_threshold_fraction <= 1:
        raise ValueError(
            "score_threshold_fraction must be between 0 and 1."
        )

    if top_n_hits is not None:
        if not isinstance(top_n_hits, int):
            raise ValueError(
                "top_n_hits must be None or an integer."
            )

        if top_n_hits < 1:
            raise ValueError(
                "top_n_hits must be at least 1."
            )

    return clean_sequence


# ============================================================
# SEQUENCE FUNCTIONS
# ============================================================


def complement_base(base):
    """
    Return the complementary DNA base.
    """

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    return complement[base]


def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA sequence.
    """

    translation = str.maketrans(
        "ACGTacgt",
        "TGCAtgca"
    )

    return sequence.translate(translation)[::-1]


# ============================================================
# PWM SCORING
# ============================================================


def score_window_with_pwm(
    window,
    pwm,
    background,
    pseudocount=1e-9
):
    """
    Score one DNA window using PWM log2 odds.

    The score at each motif position is:

        log2(PWM probability / background probability)

    Returns
    -------
    total_score : float
        Sum of all positional log-odds scores.

    positional_scores : list
        Contribution of each observed nucleotide.
    """

    positional_scores = []

    for position, base in enumerate(window):

        pwm_probability = float(
            pwm.loc[position, base]
        )

        background_probability = float(
            background[base]
        )

        score = np.log2(
            (pwm_probability + pseudocount)
            / background_probability
        )

        positional_scores.append(
            float(score)
        )

    total_score = float(
        np.sum(positional_scores)
    )

    return total_score, positional_scores


def calculate_theoretical_scores(
    pwm,
    background,
    pseudocount=1e-9
):
    """
    Calculate the maximum possible score at every PWM position.

    Also calculate the theoretical maximum score of the complete motif.
    """

    positional_maxima = []

    for position in range(len(pwm)):

        possible_scores = []

        for base in "ACGT":

            pwm_probability = float(
                pwm.loc[position, base]
            )

            background_probability = float(
                background[base]
            )

            score = np.log2(
                (pwm_probability + pseudocount)
                / background_probability
            )

            possible_scores.append(
                float(score)
            )

        positional_maxima.append(
            max(possible_scores)
        )

    theoretical_maximum = float(
        np.sum(positional_maxima)
    )

    return positional_maxima, theoretical_maximum


# ============================================================
# MOTIF SCANNING
# ============================================================


def scan_all_hits(
    sequence,
    pwm,
    background,
    threshold_fraction,
    top_n=None,
    scan_reverse=True
):
    """
    Scan all possible windows on both strands.

    Every hit with a score greater than or equal to the threshold
    is retained.
    """

    motif_length = len(pwm)

    (
        positional_maxima,
        theoretical_maximum
    ) = calculate_theoretical_scores(
        pwm,
        background
    )

    absolute_threshold = (
        threshold_fraction
        * theoretical_maximum
    )

    candidates = []

    number_of_windows = (
        len(sequence)
        - motif_length
        + 1
    )

    for start in range(number_of_windows):

        end = start + motif_length

        genomic_window = sequence[start:end]

        # ====================================================
        # POSITIVE-STRAND SCAN
        # ====================================================

        (
            plus_score,
            plus_positional_scores
        ) = score_window_with_pwm(
            genomic_window,
            pwm,
            background
        )

        candidates.append(
            {
                "start": start,
                "end": end,
                "strand": "+",
                "score": plus_score,
                "genomic_sequence": genomic_window,
                "scored_sequence": genomic_window,
                "per_position_scores": plus_positional_scores
            }
        )

        # ====================================================
        # NEGATIVE-STRAND SCAN
        # ====================================================

        if scan_reverse:

            reverse_window = reverse_complement(
                genomic_window
            )

            (
                minus_score,
                minus_positional_scores
            ) = score_window_with_pwm(
                reverse_window,
                pwm,
                background
            )

            candidates.append(
                {
                    "start": start,
                    "end": end,
                    "strand": "-",
                    "score": minus_score,

                    # Sequence as written in the positive input.
                    "genomic_sequence": genomic_window,

                    # Reverse-complement sequence scored against
                    # the PWM.
                    "scored_sequence": reverse_window,

                    # Positional scores are initially ordered in
                    # PWM/reverse-complement orientation.
                    "per_position_scores": minus_positional_scores
                }
            )

    candidates.sort(
        key=lambda hit: hit["score"],
        reverse=True
    )

    retained_hits = [
        hit
        for hit in candidates
        if hit["score"] >= absolute_threshold
    ]

    if top_n is not None:
        retained_hits = retained_hits[:top_n]

    if not retained_hits:

        best_hit = candidates[0]

        print()
        print(
            "WARNING: no motif hit passed the selected threshold."
        )
        print(
            "The best-scoring hit will be retained instead."
        )
        print(
            f"Best score: {best_hit['score']:.3f}"
        )
        print(
            f"Threshold:  {absolute_threshold:.3f}"
        )

        retained_hits = [best_hit]

    return (
        retained_hits,
        positional_maxima,
        theoretical_maximum,
        absolute_threshold
    )


# ============================================================
# LETTER HEIGHT CALCULATION
# ============================================================


def calculate_letter_height(
    raw_position_score,
    maximum_position_score,
    total_hit_score,
    theoretical_maximum
):
    """
    Calculate the plotting height of one motif nucleotide.

    The height depends on:

    1. The score of the complete motif hit.
    2. The contribution of the base at that PWM position.
    """

    if theoretical_maximum <= 0:
        return motif_min_height

    # Strength of the complete motif match.
    hit_strength = (
        total_hit_score
        / theoretical_maximum
    )

    hit_strength = float(
        np.clip(
            hit_strength,
            0.0,
            1.0
        )
    )

    # Strength of the observed base at this motif position.
    if maximum_position_score <= 0:
        position_strength = 0.0

    else:
        position_strength = (
            raw_position_score
            / maximum_position_score
        )

        position_strength = float(
            np.clip(
                position_strength,
                0.0,
                1.0
            )
        )

    height = motif_min_height + (
        hit_strength
        * position_strength
        * (
            motif_max_height
            - motif_min_height
        )
    )

    return float(height)


# ============================================================
# BUILD LOGOMAKER MATRICES
# ============================================================


def make_empty_logo_matrix(sequence_length):
    """
    Create an empty A/C/G/T Logomaker matrix.
    """

    return pd.DataFrame(
        0.0,
        index=np.arange(sequence_length),
        columns=list("ACGT")
    )


def build_logo_matrices(
    sequence,
    hits,
    positional_maxima,
    theoretical_maximum
):
    """
    Build four matrices.

    plus_baseline_matrix
        Small positive-strand bases above zero.

    minus_baseline_matrix
        Small complementary negative-strand bases below zero.

    plus_motif_matrix
        Enlarged positive-strand motif hits above zero.

    minus_motif_matrix
        Enlarged complementary negative-strand hits below zero.
    """

    sequence_length = len(sequence)

    plus_baseline_matrix = make_empty_logo_matrix(
        sequence_length
    )

    minus_baseline_matrix = make_empty_logo_matrix(
        sequence_length
    )

    plus_motif_matrix = make_empty_logo_matrix(
        sequence_length
    )

    minus_motif_matrix = make_empty_logo_matrix(
        sequence_length
    )

    # ========================================================
    # COMPLETE POSITIVE AND NEGATIVE STRANDS
    # ========================================================

    for genomic_position, plus_base in enumerate(sequence):

        minus_base = complement_base(
            plus_base
        )

        # Positive input sequence.
        plus_baseline_matrix.loc[
            genomic_position,
            plus_base
        ] = baseline_height

        # Complementary negative strand.
        minus_baseline_matrix.loc[
            genomic_position,
            minus_base
        ] = -baseline_height

    # ========================================================
    # RETAINED MOTIF HITS
    # ========================================================

    for hit in hits:

        positional_scores = list(
            hit["per_position_scores"]
        )

        local_positional_maxima = list(
            positional_maxima
        )

        # Negative-strand scores were calculated in reverse-
        # complement orientation. Reverse the score arrays so
        # that they align with left-to-right genomic positions.
        if hit["strand"] == "-":

            positional_scores = positional_scores[::-1]

            local_positional_maxima = (
                local_positional_maxima[::-1]
            )

        for offset, raw_position_score in enumerate(
            positional_scores
        ):

            genomic_position = (
                hit["start"]
                + offset
            )

            plus_base = sequence[
                genomic_position
            ]

            minus_base = complement_base(
                plus_base
            )

            letter_height = calculate_letter_height(
                raw_position_score=raw_position_score,
                maximum_position_score=(
                    local_positional_maxima[offset]
                ),
                total_hit_score=hit["score"],
                theoretical_maximum=theoretical_maximum
            )

            # =================================================
            # POSITIVE-STRAND MOTIF HIT
            # =================================================

            if hit["strand"] == "+":

                current_height = (
                    plus_motif_matrix.loc[
                        genomic_position,
                        plus_base
                    ]
                )

                # When hits overlap, keep the tallest value.
                if letter_height > current_height:

                    plus_motif_matrix.loc[
                        genomic_position,
                        :
                    ] = 0.0

                    plus_motif_matrix.loc[
                        genomic_position,
                        plus_base
                    ] = letter_height

            # =================================================
            # NEGATIVE-STRAND MOTIF HIT
            # =================================================

            else:

                current_height = abs(
                    minus_motif_matrix.loc[
                        genomic_position,
                        minus_base
                    ]
                )

                # Negative values place letters below zero.
                if letter_height > current_height:

                    minus_motif_matrix.loc[
                        genomic_position,
                        :
                    ] = 0.0

                    minus_motif_matrix.loc[
                        genomic_position,
                        minus_base
                    ] = -letter_height

    return (
        plus_baseline_matrix,
        minus_baseline_matrix,
        plus_motif_matrix,
        minus_motif_matrix
    )


# ============================================================
# TERMINAL OUTPUT
# ============================================================


def print_hits(
    hits,
    absolute_threshold,
    theoretical_maximum
):
    """
    Print all retained hits in the Terminal.
    """

    print()
    print("PWM motif scan")
    print(
        "======================================================================"
    )
    print(
        f"Threshold fraction:       "
        f"{score_threshold_fraction:.3f}"
    )
    print(
        f"Absolute score threshold: "
        f"{absolute_threshold:.3f}"
    )
    print(
        f"Theoretical maximum:      "
        f"{theoretical_maximum:.3f}"
    )
    print(
        f"Retained hits:            "
        f"{len(hits)}"
    )
    print(
        "======================================================================"
    )
    print(
        "N\tstart\tend\tstrand\tscore\t"
        "genomic\tPWM_orientation"
    )

    for hit_number, hit in enumerate(
        hits,
        start=1
    ):

        print(
            f"{hit_number}\t"
            f"{hit['start'] + 1}\t"
            f"{hit['end']}\t"
            f"{hit['strand']}\t"
            f"{hit['score']:.3f}\t"
            f"{hit['genomic_sequence']}\t"
            f"{hit['scored_sequence']}"
        )

    print(
        "======================================================================"
    )


# ============================================================
# PLOTTING
# ============================================================


def plot_logo(
    sequence,
    hits,
    plus_baseline_matrix,
    minus_baseline_matrix,
    plus_motif_matrix,
    minus_motif_matrix,
    absolute_threshold,
    theoretical_maximum
):
    """
    Plot the complete strand-aware motif figure.
    """

    figure_width = max(
        minimum_figure_width,
        len(sequence) * width_per_base
    )

    fig, ax = plt.subplots(
        figsize=(
            figure_width,
            figure_height
        )
    )

    # ========================================================
    # HIGHLIGHT RETAINED MOTIF REGIONS
    # ========================================================

    for hit in hits:

        if hit["strand"] == "+":
            highlight_color = plus_highlight_color

        else:
            highlight_color = minus_highlight_color

        ax.axvspan(
            hit["start"] - 0.5,
            hit["end"] - 0.5,
            color=highlight_color,
            alpha=highlight_alpha,
            zorder=0
        )

    # ========================================================
    # COMPLETE POSITIVE-STRAND BASELINE
    # ========================================================

    logomaker.Logo(
        plus_baseline_matrix,
        ax=ax,
        color_scheme="classic",
        stack_order="big_on_top"
    )

    # ========================================================
    # COMPLETE NEGATIVE-STRAND BASELINE
    # ========================================================

    logomaker.Logo(
        minus_baseline_matrix,
        ax=ax,
        color_scheme="classic",
        stack_order="big_on_top",
        flip_below=True
    )

    # ========================================================
    # POSITIVE-STRAND MOTIF HITS
    # ========================================================

    logomaker.Logo(
        plus_motif_matrix,
        ax=ax,
        color_scheme="classic",
        stack_order="big_on_top"
    )

    # ========================================================
    # NEGATIVE-STRAND MOTIF HITS
    # ========================================================

    logomaker.Logo(
        minus_motif_matrix,
        ax=ax,
        color_scheme="classic",
        stack_order="big_on_top",
        flip_below=True
    )

    # Central zero line.
    ax.axhline(
        y=0,
        linewidth=1.0,
        color="black",
        zorder=20
    )

    # ========================================================
    # COUNTS AND LABELS
    # ========================================================

    plus_hit_count = sum(
        hit["strand"] == "+"
        for hit in hits
    )

    minus_hit_count = sum(
        hit["strand"] == "-"
        for hit in hits
    )

    title = (
        f"PWM motif scan: "
        f"{plus_hit_count} positive-strand hit(s), "
        f"{minus_hit_count} negative-strand hit(s)"
    )

    ax.set_title(
        title,
        fontsize=13,
        pad=24
    )

    subtitle = (
        f"Threshold = "
        f"{score_threshold_fraction:.2f} × theoretical maximum "
        f"({absolute_threshold:.2f}/{theoretical_maximum:.2f})"
    )

    ax.text(
        0.01,
        1.015,
        subtitle,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="bottom"
    )

    ax.set_ylabel(
        "PWM-weighted letter height"
    )

    ax.set_xlabel(
        "Position in input DNA sequence"
    )

    sequence_positions = np.arange(
        len(sequence)
    )

    ax.set_xticks(
        sequence_positions
    )

    ax.set_xticklabels(
        sequence_positions + 1,
        rotation=90,
        fontsize=7
    )

    ax.set_xlim(
        -0.5,
        len(sequence) - 0.5
    )

    # Symmetrical y-axis.
    y_axis_limit = (
        motif_max_height
        * 1.12
    )

    ax.set_ylim(
        -y_axis_limit,
        y_axis_limit
    )

    # Strand labels.
    ax.text(
        -0.015,
        0.76,
        "+",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        horizontalalignment="right",
        verticalalignment="center"
    )

    ax.text(
        -0.015,
        0.24,
        "−",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        horizontalalignment="right",
        verticalalignment="center"
    )

    ax.spines["top"].set_visible(
        False
    )

    ax.spines["right"].set_visible(
        False
    )

    plt.tight_layout()

    # ========================================================
    # EXPORT
    # ========================================================

    fig.savefig(
        svg_output,
        format="svg",
        bbox_inches="tight"
    )

    fig.savefig(
        pdf_output,
        format="pdf",
        bbox_inches="tight"
    )

    fig.savefig(
        png_output,
        format="png",
        dpi=300,
        bbox_inches="tight"
    )

    print()
    print("Output files created:")
    print(
        f"SVG: {os.path.abspath(svg_output)}"
    )
    print(
        f"PDF: {os.path.abspath(pdf_output)}"
    )
    print(
        f"PNG: {os.path.abspath(png_output)}"
    )
    print()

    if show_plot:
        plt.show()

    else:
        plt.close(fig)


# ============================================================
# MAIN
# ============================================================


def main():
    """
    Run the complete analysis and plotting workflow.
    """

    clean_sequence = clean_and_validate_inputs(
        sequence,
        pwm
    )

    (
        hits,
        positional_maxima,
        theoretical_maximum,
        absolute_threshold
    ) = scan_all_hits(
        sequence=clean_sequence,
        pwm=pwm,
        background=background,
        threshold_fraction=score_threshold_fraction,
        top_n=top_n_hits,
        scan_reverse=scan_reverse_strand
    )

    print_hits(
        hits=hits,
        absolute_threshold=absolute_threshold,
        theoretical_maximum=theoretical_maximum
    )

    (
        plus_baseline_matrix,
        minus_baseline_matrix,
        plus_motif_matrix,
        minus_motif_matrix
    ) = build_logo_matrices(
        sequence=clean_sequence,
        hits=hits,
        positional_maxima=positional_maxima,
        theoretical_maximum=theoretical_maximum
    )

    plot_logo(
        sequence=clean_sequence,
        hits=hits,
        plus_baseline_matrix=plus_baseline_matrix,
        minus_baseline_matrix=minus_baseline_matrix,
        plus_motif_matrix=plus_motif_matrix,
        minus_motif_matrix=minus_motif_matrix,
        absolute_threshold=absolute_threshold,
        theoretical_maximum=theoretical_maximum
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print()
        print("ERROR:")
        print(error)
        print()

        sys.exit(1)