#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def sample_name_from_count_column(colname: str) -> str:
    name = Path(colname).name
    for suffix in [
        ".Aligned.sortedByCoord.out.bam",
        ".sortedByCoord.out.bam",
        ".bam",
    ]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_featurecounts_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")

    if "Geneid" not in df.columns:
        raise ValueError(f"{path} does not contain a 'Geneid' column")

    count_col = df.columns[-1]
    sample = sample_name_from_count_column(count_col)

    out = df[["Geneid", count_col]].copy()
    out = out.rename(columns={count_col: sample})
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge per-sample featureCounts TSV files into one gene count matrix."
    )
    parser.add_argument(
        "--counts",
        nargs="+",
        required=True,
        help="Per-sample featureCounts .tsv files.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output merged count matrix TSV.",
    )
    args = parser.parse_args()

    count_files = [Path(p) for p in args.counts]

    if not count_files:
        raise ValueError("No featureCounts files provided")

    merged = None

    for count_file in count_files:
        sample_df = read_featurecounts_file(count_file)
        if merged is None:
            merged = sample_df
        else:
            merged = merged.merge(sample_df, on="Geneid", how="outer")

    merged = merged.fillna(0)

    sample_cols = [c for c in merged.columns if c != "Geneid"]
    merged[sample_cols] = merged[sample_cols].astype(int)

    merged.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()