#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path


# Match decimal numbers, including scientific notation.
NUMBER = r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"

# Paired-end: read 1 agrees with transcript orientation.
PE_FORWARD_PATTERN = re.compile(
    r'Fraction of reads explained by "1\+\+,1--,2\+\-,2-\+":\s*'
    + NUMBER
)

# Paired-end: read 1 opposes transcript orientation.
PE_REVERSE_PATTERN = re.compile(
    r'Fraction of reads explained by "1\+\-,1-\+,2\+\+,2--":\s*'
    + NUMBER
)

# Single-end: read agrees with transcript orientation.
SE_FORWARD_PATTERN = re.compile(
    r'Fraction of reads explained by "\+\+,--":\s*'
    + NUMBER
)

# Single-end: read opposes transcript orientation.
SE_REVERSE_PATTERN = re.compile(
    r'Fraction of reads explained by "\+-,\s*-\+":\s*'
    + NUMBER
)


def extract_fractions(text: str) -> tuple[float, float]:

    # Try paired-end format first
    forward_match = PE_FORWARD_PATTERN.search(text)
    reverse_match = PE_REVERSE_PATTERN.search(text)

    if forward_match and reverse_match:
        return (
            float(forward_match.group(1)),
            float(reverse_match.group(1)),
        )

    # Try single-end format
    forward_match = SE_FORWARD_PATTERN.search(text)
    reverse_match = SE_REVERSE_PATTERN.search(text)

    if forward_match and reverse_match:
        return (
            float(forward_match.group(1)),
            float(reverse_match.group(1)),
        )

    raise ValueError(
        "Could not identify strandedness fractions "
        "in the RSeQC report."
    )


def classify_report(
    path: Path,
    threshold: float,
    unstranded_tolerance: float,
) -> str:

    text = path.read_text(encoding="utf-8")

    forward_fraction, reverse_fraction = extract_fractions(text)

    if not (0 <= forward_fraction <= 1):
        raise ValueError("Invalid forward fraction.")

    if not (0 <= reverse_fraction <= 1):
        raise ValueError("Invalid reverse fraction.")

    total = forward_fraction + reverse_fraction

    if total <= 0:
        raise ValueError(
            "No informative reads available for strandedness inference."
        )

    # Relative proportions among informative reads
    forward_relative = forward_fraction / total
    reverse_relative = reverse_fraction / total

    if (
        reverse_relative >= threshold
        and reverse_relative > forward_relative
    ):
        return "REVERSE"

    if (
        forward_relative >= threshold
        and forward_relative > reverse_relative
    ):
        return "FORWARD"

    # Approximately equal orientations indicate unstranded RNA-seq
    if abs(forward_relative - 0.5) <= unstranded_tolerance:
        return "UNSTRANDED"

    raise ValueError(
        f"Ambiguous strandedness in {path.name}: "
        f"forward={forward_relative:.3f}, "
        f"reverse={reverse_relative:.3f}. "
        "Neither orientation reaches the confidence threshold, "
        "and the library is not sufficiently balanced to classify "
        "as unstranded."
    )


def main() -> None:

    parser = argparse.ArgumentParser(
        description="Infer RNA-seq library strandedness from one RSeQC report."
    )

    parser.add_argument(
        "report",
        type=Path,
        help="RSeQC infer_experiment.py output file.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.8,
        help="Minimum relative orientation fraction for confident strandedness.",
    )

    parser.add_argument(
        "--unstranded-tolerance",
        type=float,
        default=0.15,
        help="Maximum deviation from 0.5 for an unstranded classification.",
    )

    args = parser.parse_args()

    if not 0.5 < args.threshold <= 1:
        parser.error("--threshold must be greater than 0.5 and at most 1.")

    if not 0 <= args.unstranded_tolerance < 0.5:
        parser.error(
            "--unstranded-tolerance must be between 0 and 0.5."
        )

    result = classify_report(
        args.report,
        args.threshold,
        args.unstranded_tolerance,
    )

    print(result)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)