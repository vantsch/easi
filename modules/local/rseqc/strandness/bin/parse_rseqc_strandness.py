#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path


FORWARD_PATTERN = re.compile(
    r'Fraction of reads explained by "1\+\+,1--,2\+\-,2-\+":\s*([0-9.]+)'
)

REVERSE_PATTERN = re.compile(
    r'Fraction of reads explained by "1\+\-,1-\+,2\+\+,2--":\s*([0-9.]+)'
)


def classify_report(path: Path, threshold: float) -> str:
    text = path.read_text(encoding="utf-8")

    forward_match = FORWARD_PATTERN.search(text)
    reverse_match = REVERSE_PATTERN.search(text)

    if not forward_match or not reverse_match:
        raise ValueError(
            f"Could not find paired-end strandness fractions in {path}"
        )

    forward_fraction = float(forward_match.group(1))
    reverse_fraction = float(reverse_match.group(1))

    if reverse_fraction >= threshold and reverse_fraction > forward_fraction:
        return "REVERSE"

    if forward_fraction >= threshold and forward_fraction > reverse_fraction:
        return "FORWARD"

    return "UNSTRANDED"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reports", nargs="+", type=Path)
    parser.add_argument("--threshold", type=float, default=0.8)
    args = parser.parse_args()

    classifications = {
        report.name: classify_report(report, args.threshold)
        for report in args.reports
    }

    unique = set(classifications.values())

    if len(unique) != 1:
        details = ", ".join(
            f"{sample}={classification}"
            for sample, classification in classifications.items()
        )
        raise RuntimeError(
            f"RSeQC strandness classifications disagree across samples: {details}"
        )

    print(unique.pop())


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)