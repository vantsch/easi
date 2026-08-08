#!/usr/bin/env python3

"""
Integrate significant splicing regions reported by MAJIQ, rMATS and LeafCutter.

The primary comparison unit is a tool-level event region:

  MAJIQ       module_id
  rMATS       event_id
  LeafCutter  cluster

A region spans the minimum intron_start to maximum intron_end of all
connections assigned to that event. This gives the three tools a common
region definition based on splice connections rather than tool-specific exon
footprints.

A significant multi-tool overlap requires:

  1. same comparison scope, chromosome and strand;
  2. both tool events are significant according to their native statistic;
  3. reciprocal regional overlap >= --overlap-threshold (default 0.50).

Exact connections, retained-intron/excised-intron matches and shared donor or
acceptor sites are reported as supporting annotations. They are not required
for a regional hit.

Outputs
-------
<prefix>.normalized_connections.tsv
<prefix>.tool_event_regions.tsv
<prefix>.all_region_overlaps.tsv
<prefix>.significant_region_overlaps.tsv
<prefix>.multi_tool_region_summary.tsv
<prefix>.ranked_splicing_events.tsv
<prefix>.splice_region_integration.xlsx
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


CANONICAL_COLUMNS = [
    "dataset_id",
    "contrast",
    "tool",
    "tool_event_id",
    "tool_group_id",
    "connection_role",
    "gene_id",
    "gene_name",
    "chromosome",
    "strand",
    "intron_start",
    "intron_end",
    "feature_class",
    "event_type",
    "effect_value",
    "effect_metric",
    "event_significant",
    "connection_significant",
    "significance_value",
    "significance_metric",
    "source_file",
    "source_row",
]


SUPPORT_PRIORITY = {
    "EXACT_CONNECTION": 1,
    "SAME_INTRON_LOCUS": 2,
    "SHARED_SPLICE_SITE": 3,
    "REGIONAL_OVERLAP": 4,
}


EXCEL_MAX_DATA_ROWS = 1_048_575


class UnionFind:

    def __init__(
        self,
        values: Iterable[str],
    ) -> None:
        self.parent = {
            value: value
            for value in values
        }

        self.rank = {
            value: 0
            for value in values
        }

    def find(
        self,
        value: str,
    ) -> str:
        if self.parent[value] != value:
            self.parent[value] = self.find(
                self.parent[value]
            )

        return self.parent[value]

    def union(
        self,
        left: str,
        right: str,
    ) -> None:
        left_root = self.find(left)
        right_root = self.find(right)

        if left_root == right_root:
            return

        if (
            self.rank[left_root]
            < self.rank[right_root]
        ):
            left_root, right_root = (
                right_root,
                left_root,
            )

        self.parent[right_root] = left_root

        if (
            self.rank[left_root]
            == self.rank[right_root]
        ):
            self.rank[left_root] += 1


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Identify significant, reciprocally overlapping "
            "splicing regions reported by MAJIQ, rMATS "
            "and LeafCutter."
        )
    )

    parser.add_argument(
        "--majiq-junctions",
        required=True,
        type=Path,
        help="MAJIQ delta-PSI modulize junctions.tsv.",
    )

    parser.add_argument(
        "--rmats-jc",
        required=True,
        nargs="+",
        type=Path,
        help="rMATS *.MATS.JC.annotated.tsv files.",
    )

    parser.add_argument(
        "--leafcutter-introns",
        required=True,
        type=Path,
        help="LeafCutter *.leafcutter.introns.annotated.tsv.",
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Output prefix, normally the contrast name.",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Output directory.",
    )

    parser.add_argument(
        "--overlap-threshold",
        type=float,
        default=0.50,
        help=(
            "Minimum fraction of both event regions that "
            "must overlap. Default: 0.50."
        ),
    )

    parser.add_argument(
        "--sashimi-top-n",
        type=int,
        default=30,
        help=(
            "Number of highest-ranked significant regions to mark "
            "for sashimi plotting. Default: 30."
        ),
    )

    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help=(
            "FDR/adjusted-P threshold for rMATS and "
            "LeafCutter."
        ),
    )

    parser.add_argument(
        "--majiq-probability-threshold",
        type=float,
        default=0.95,
        help=(
            "MAJIQ posterior probability threshold when "
            "no explicit flag is present."
        ),
    )

    args = parser.parse_args()

    if not 0 < args.overlap_threshold <= 1:
        parser.error(
            "--overlap-threshold must be > 0 and <= 1"
        )

    if args.sashimi_top_n < 1:
        parser.error(
            "--sashimi-top-n must be at least 1"
        )

    return args


def normalise_column_name(
    value: str,
) -> str:
    value = str(value).strip().lower()

    value = re.sub(
        r"[^a-z0-9]+",
        "_",
        value,
    )

    return value.strip("_")


def read_tsv(
    path: Path,
) -> pd.DataFrame:
    """
    Read a TSV while tolerating metadata lines before the real header.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"Input file does not exist: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        lines = handle.readlines()

    if not lines:
        raise ValueError(
            f"Input file is empty: {path}"
        )

    header_tokens = {
        "module_id",
        "gene_id",
        "geneid",
        "event_id",
        "intron_id",
        "junction_coord",
        "cluster",
        "seqname",
        "seqid",
        "chr",
        "strand",
        "fdr",
        "contrast",
    }

    header_index: int | None = None
    fallback_index: int | None = None

    for index, line in enumerate(lines):
        stripped = line.rstrip(
            "\r\n"
        )

        if not stripped:
            continue

        fields = stripped.split(
            "\t"
        )

        if len(fields) < 2:
            continue

        normalised_fields = {
            normalise_column_name(field)
            for field in fields
        }

        if (
            len(
                normalised_fields
                & header_tokens
            )
            >= 2
        ):
            header_index = index
            break

        if fallback_index is None:
            for next_index in range(
                index + 1,
                len(lines),
            ):
                next_line = lines[
                    next_index
                ].rstrip(
                    "\r\n"
                )

                if not next_line:
                    continue

                if (
                    len(
                        next_line.split(
                            "\t"
                        )
                    )
                    == len(fields)
                ):
                    fallback_index = index

                break

    if header_index is None:
        header_index = fallback_index

    if header_index is None:
        raise ValueError(
            f"Could not identify a tabular header in {path}"
        )

    frame = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
        skiprows=header_index,
    )

    frame.columns = [
        normalise_column_name(column)
        for column in frame.columns
    ]

    if frame.columns.duplicated().any():
        duplicates = sorted(
            set(
                frame.columns[
                    frame.columns.duplicated()
                ]
            )
        )

        raise ValueError(
            f"{path} has duplicate columns after normalisation: "
            f"{', '.join(duplicates)}"
        )

    frame["source_file"] = str(path)

    frame["source_row"] = range(
        header_index + 2,
        header_index + 2 + len(frame),
    )

    print(
        f"Reading {path}: "
        f"header at line {header_index + 1}, "
        f"{len(frame.columns) - 2} columns"
    )

    return frame


def find_column(
    frame: pd.DataFrame,
    candidates: Iterable[str],
    *,
    required: bool = False,
    label: str = "column",
    path: Path | None = None,
) -> str | None:

    for candidate in candidates:
        normalised = normalise_column_name(
            candidate
        )

        if normalised in frame.columns:
            return normalised

    if required:
        location = (
            f" in {path}"
            if path is not None
            else ""
        )

        raise ValueError(
            f"Could not identify {label}{location}. "
            f"Available columns: "
            f"{', '.join(frame.columns)}"
        )

    return None


def find_matching_column(
    frame: pd.DataFrame,
    predicate: Any,
) -> str | None:

    return next(
        (
            column
            for column in frame.columns
            if predicate(column)
        ),
        None,
    )


def require_columns(
    frame: pd.DataFrame,
    path: Path,
    columns: Iterable[str],
) -> None:

    missing = [
        normalise_column_name(column)
        for column in columns
        if normalise_column_name(column)
        not in frame.columns
    ]

    if missing:
        raise ValueError(
            f"{path} is missing columns: "
            f"{', '.join(missing)}. "
            f"Available columns: "
            f"{', '.join(frame.columns)}"
        )


def clean_text(
    series: pd.Series,
) -> pd.Series:

    return (
        series.astype("string")
        .str.strip()
        .replace(
            {
                "": pd.NA,
                "NA": pd.NA,
                "N/A": pd.NA,
                "nan": pd.NA,
                "None": pd.NA,
                "null": pd.NA,
            }
        )
    )


def parse_boolean(
    value: Any,
) -> Any:

    if pd.isna(value):
        return pd.NA

    text = normalise_column_name(
        str(value)
    )

    if text in {
        "true",
        "t",
        "yes",
        "y",
        "1",
        "significant",
        "sig",
        "changing",
        "tested_significant",
    }:
        return True

    if text in {
        "false",
        "f",
        "no",
        "n",
        "0",
        "not_significant",
        "nonsignificant",
        "ns",
        "non_changing",
    }:
        return False

    return pd.NA


def boolean_series(
    frame: pd.DataFrame,
    column: str | None,
) -> pd.Series:

    if column is None:
        return pd.Series(
            pd.NA,
            index=frame.index,
            dtype="boolean",
        )

    return (
        frame[column]
        .map(parse_boolean)
        .astype("boolean")
    )


def numeric_series(
    frame: pd.DataFrame,
    column: str | None,
) -> pd.Series:

    if column is None:
        return pd.Series(
            pd.NA,
            index=frame.index,
            dtype="Float64",
        )

    return pd.to_numeric(
        frame[column],
        errors="coerce",
    )


def fill_significance(
    explicit: pd.Series,
    statistic: pd.Series,
    *,
    mode: str,
    threshold: float,
) -> pd.Series:

    result = explicit.astype(
        "boolean"
    ).copy()

    missing = result.isna()

    if mode == "maximum":
        result.loc[missing] = (
            statistic.loc[missing]
            >= threshold
        )

    elif mode == "minimum":
        result.loc[missing] = (
            statistic.loc[missing]
            <= threshold
        )

    else:
        raise ValueError(
            f"Unsupported significance mode: {mode}"
        )

    result.loc[
        statistic.isna()
        & explicit.isna()
    ] = pd.NA

    return result.astype(
        "boolean"
    )


def parse_junction_coordinates(
    series: pd.Series,
) -> pd.DataFrame:

    patterns = [
        re.compile(
            r"(?:(?P<chrom>[^:\s]+):)?"
            r"(?P<start>\d[\d,]*)"
            r"[-:]"
            r"(?P<end>\d[\d,]*)"
        ),
        re.compile(
            r"(?:(?P<chrom>[^_\s]+)_)?"
            r"(?P<start>\d[\d,]*)"
            r"_"
            r"(?P<end>\d[\d,]*)"
        ),
    ]

    rows: list[dict[str, Any]] = []

    for value in series:
        parsed: dict[str, Any] = {
            "chromosome": pd.NA,
            "intron_start": pd.NA,
            "intron_end": pd.NA,
        }

        if pd.notna(value):
            text = str(value).strip()

            for pattern in patterns:
                match = pattern.search(
                    text
                )

                if match is None:
                    continue

                parsed = {
                    "chromosome": (
                        match.groupdict()
                        .get(
                            "chrom",
                            pd.NA,
                        )
                    ),
                    "intron_start": int(
                        match.group(
                            "start"
                        ).replace(
                            ",",
                            "",
                        )
                    ),
                    "intron_end": int(
                        match.group(
                            "end"
                        ).replace(
                            ",",
                            "",
                        )
                    ),
                }

                break

        rows.append(
            parsed
        )

    return pd.DataFrame(
        rows,
        index=series.index,
    )


def canonical_frame(
    rows: list[dict[str, Any]],
) -> pd.DataFrame:

    frame = pd.DataFrame(
        rows
    )

    for column in CANONICAL_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA

    return frame[
        CANONICAL_COLUMNS
    ]


def parse_majiq(
    path: Path,
    contrast: str,
    *,
    probability_threshold: float,
) -> pd.DataFrame:

    frame = read_tsv(
        path
    )

    module_column = find_column(
        frame,
        [
            "module_id",
            "module",
        ],
        required=True,
        label="MAJIQ module ID",
        path=path,
    )

    event_column = find_column(
        frame,
        [
            "event_id",
            "lsv_id",
            "junction_name",
        ],
        required=True,
        label="MAJIQ event ID",
        path=path,
    )

    coordinate_column = find_column(
        frame,
        [
            "junction_coord",
            "junction_coordinate",
            "coordinates",
        ],
        required=True,
        label="MAJIQ junction coordinates",
        path=path,
    )

    chromosome_column = find_column(
        frame,
        [
            "seqid",
            "seqname",
            "chromosome",
            "chrom",
            "chr",
        ],
        required=True,
        label="MAJIQ chromosome",
        path=path,
    )

    strand_column = find_column(
        frame,
        [
            "strand",
            "gene_strand",
            "lsv_strand",
        ],
        required=True,
        label="MAJIQ strand",
        path=path,
    )

    parsed_coordinates = (
        parse_junction_coordinates(
            frame[coordinate_column]
        )
    )

    starts = pd.to_numeric(
        parsed_coordinates[
            "intron_start"
        ],
        errors="coerce",
    )

    ends = pd.to_numeric(
        parsed_coordinates[
            "intron_end"
        ],
        errors="coerce",
    )

    gene_id_column = find_column(
        frame,
        [
            "gene_id",
            "geneid",
        ],
    )

    gene_name_column = find_column(
        frame,
        [
            "gene_name",
            "gene_symbol",
            "genesymbol",
            "symbol",
        ],
    )

    event_type_column = find_column(
        frame,
        [
            "module_event_combination",
            "event_type",
            "module_type",
            "lsv_type",
        ],
    )

    junction_name_column = find_column(
        frame,
        [
            "junction_name",
            "junction_id",
        ],
    )

    contrast_column = find_column(
        frame,
        [
            "contrast",
            "comparison",
        ],
    )

    dataset_column = find_column(
        frame,
        [
            "dataset_id",
            "dataset",
        ],
    )

    effect_column = find_column(
        frame,
        [
            "dpsi",
            "delta_psi",
            "median_dpsi",
            "mean_dpsi",
        ],
    )

    if effect_column is None:
        effect_column = find_matching_column(
            frame,
            lambda column: (
                column.endswith(
                    "_median_dpsi"
                )
                or column.endswith(
                    "_mean_dpsi"
                )
                or column.endswith(
                    "_dpsi"
                )
            ),
        )

    probability_column = find_column(
        frame,
        [
            "probability_changing",
            "posterior_probability",
        ],
    )

    if probability_column is None:
        probability_column = find_matching_column(
            frame,
            lambda column: (
                column.endswith(
                    "_probability_changing"
                )
                or column.endswith(
                    "_posterior_probability"
                )
            ),
        )

    event_flag_column = find_column(
        frame,
        [
            "event_changing",
            "significant",
            "is_significant",
        ],
    )

    junction_flag_column = find_column(
        frame,
        [
            "junction_changing",
            "junction_significant",
        ],
    )

    probabilities = numeric_series(
        frame,
        probability_column,
    )

    event_significant = fill_significance(
        boolean_series(
            frame,
            event_flag_column,
        ),
        probabilities,
        mode="maximum",
        threshold=probability_threshold,
    )

    connection_significant = boolean_series(
        frame,
        junction_flag_column,
    )

    connection_significant = (
        connection_significant.fillna(
            event_significant
        )
    )

    rows: list[dict[str, Any]] = []

    for index in frame.index:
        rows.append(
            {
                "dataset_id": (
                    frame.at[
                        index,
                        dataset_column,
                    ]
                    if dataset_column is not None
                    else pd.NA
                ),
                "contrast": (
                    frame.at[
                        index,
                        contrast_column,
                    ]
                    if contrast_column is not None
                    else contrast
                ),
                "tool": "MAJIQ",
                "tool_event_id": (
                    frame.at[
                        index,
                        event_column,
                    ]
                ),
                "tool_group_id": (
                    frame.at[
                        index,
                        module_column,
                    ]
                ),
                "connection_role": (
                    frame.at[
                        index,
                        junction_name_column,
                    ]
                    if junction_name_column is not None
                    else "junction"
                ),
                "gene_id": (
                    frame.at[
                        index,
                        gene_id_column,
                    ]
                    if gene_id_column is not None
                    else pd.NA
                ),
                "gene_name": (
                    frame.at[
                        index,
                        gene_name_column,
                    ]
                    if gene_name_column is not None
                    else pd.NA
                ),
                "chromosome": (
                    frame.at[
                        index,
                        chromosome_column,
                    ]
                ),
                "strand": (
                    frame.at[
                        index,
                        strand_column,
                    ]
                ),
                "intron_start": (
                    starts.at[index]
                ),
                "intron_end": (
                    ends.at[index]
                ),
                "feature_class": (
                    "EXCISED_INTRON"
                ),
                "event_type": (
                    frame.at[
                        index,
                        event_type_column,
                    ]
                    if event_type_column is not None
                    else "LSV"
                ),
                "effect_value": (
                    frame.at[
                        index,
                        effect_column,
                    ]
                    if effect_column is not None
                    else pd.NA
                ),
                "effect_metric": (
                    effect_column
                    or pd.NA
                ),
                "event_significant": (
                    event_significant.at[
                        index
                    ]
                ),
                "connection_significant": (
                    connection_significant.at[
                        index
                    ]
                ),
                "significance_value": (
                    probabilities.at[
                        index
                    ]
                ),
                "significance_metric": (
                    probability_column
                    or pd.NA
                ),
                "source_file": (
                    frame.at[
                        index,
                        "source_file",
                    ]
                ),
                "source_row": (
                    frame.at[
                        index,
                        "source_row",
                    ]
                ),
            }
        )

    return canonical_frame(
        rows
    )


def infer_rmats_event_type(
    path: Path,
) -> str:

    upper_name = path.name.upper()

    for event_type in [
        "SE",
        "A3SS",
        "A5SS",
        "MXE",
        "RI",
    ]:
        if re.search(
            rf"(?:^|[._-])"
            rf"{event_type}"
            rf"(?:[._-]|$)",
            upper_name,
        ):
            return event_type

    raise ValueError(
        f"Could not infer rMATS event type "
        f"from {path.name}"
    )


def format_rmats_event_id(
    event_type: str,
    value: Any,
) -> str:

    text = str(value)

    if text.upper().startswith(
        f"{event_type}:"
    ):
        return text

    return f"{event_type}:{text}"


def parse_rmats_file(
    path: Path,
    contrast: str,
    *,
    fdr_threshold: float,
) -> pd.DataFrame:

    frame = read_tsv(
        path
    )

    event_type = infer_rmats_event_type(
        path
    )

    required_event_columns = {
        "SE": [
            "exonStart_0base",
            "exonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE",
        ],
        "A3SS": [
            "longExonStart_0base",
            "longExonEnd",
            "shortES",
            "shortEE",
            "flankingES",
            "flankingEE",
        ],
        "A5SS": [
            "longExonStart_0base",
            "longExonEnd",
            "shortES",
            "shortEE",
            "flankingES",
            "flankingEE",
        ],
        "MXE": [
            "1stExonStart_0base",
            "1stExonEnd",
            "2ndExonStart_0base",
            "2ndExonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE",
        ],
        "RI": [
            "riExonStart_0base",
            "riExonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE",
        ],
    }

    event_id_column = find_column(
        frame,
        [
            "event_id",
            "id",
            "event_id_2",
        ],
        required=True,
        label="rMATS event ID",
        path=path,
    )

    effect_column = find_column(
        frame,
        [
            "dpsi_target_minus_reference",
            "incleveldifference",
        ],
        required=True,
        label="rMATS effect size",
        path=path,
    )

    fdr_column = find_column(
        frame,
        [
            "fdr",
            "q_value",
            "qvalue",
            "padj",
        ],
        required=True,
        label="rMATS FDR",
        path=path,
    )

    require_columns(
        frame,
        path,
        [
            "chr",
            "strand",
        ]
        + required_event_columns[
            event_type
        ],
    )

    gene_id_column = find_column(
        frame,
        [
            "geneid",
            "gene_id",
        ],
    )

    gene_name_column = find_column(
        frame,
        [
            "genesymbol",
            "gene_symbol",
            "gene_name",
            "symbol",
        ],
    )

    contrast_column = find_column(
        frame,
        [
            "contrast",
            "comparison",
        ],
    )

    dataset_column = find_column(
        frame,
        [
            "dataset_id",
            "dataset",
        ],
    )

    def coordinate(
        row: pd.Series,
        column: str,
    ) -> int:

        normalised = (
            normalise_column_name(
                column
            )
        )

        value = pd.to_numeric(
            row[normalised],
            errors="coerce",
        )

        if pd.isna(value):
            raise ValueError(
                f"Non-numeric coordinate in {path}, "
                f"row {row['source_row']}, "
                f"column {normalised}"
            )

        return int(value)

    rows: list[dict[str, Any]] = []

    for _, row in frame.iterrows():
        strand = str(
            row["strand"]
        ).strip()

        event_id = format_rmats_event_id(
            event_type,
            row[event_id_column],
        )

        effect = pd.to_numeric(
            row[effect_column],
            errors="coerce",
        )

        fdr = pd.to_numeric(
            row[fdr_column],
            errors="coerce",
        )

        significant = (
            bool(
                fdr <= fdr_threshold
            )
            if pd.notna(fdr)
            else pd.NA
        )

        if event_type == "SE":
            exon_start = coordinate(
                row,
                "exonStart_0base",
            )

            exon_end = coordinate(
                row,
                "exonEnd",
            )

            upstream_ee = coordinate(
                row,
                "upstreamEE",
            )

            downstream_es = coordinate(
                row,
                "downstreamES",
            )

            connections = [
                (
                    upstream_ee,
                    exon_start,
                    "inclusion_upstream_junction",
                    "EXCISED_INTRON",
                ),
                (
                    exon_end,
                    downstream_es,
                    "inclusion_downstream_junction",
                    "EXCISED_INTRON",
                ),
                (
                    upstream_ee,
                    downstream_es,
                    "skipping_junction",
                    "EXCISED_INTRON",
                ),
            ]

        elif event_type == "MXE":
            first_start = coordinate(
                row,
                "1stExonStart_0base",
            )

            first_end = coordinate(
                row,
                "1stExonEnd",
            )

            second_start = coordinate(
                row,
                "2ndExonStart_0base",
            )

            second_end = coordinate(
                row,
                "2ndExonEnd",
            )

            upstream_ee = coordinate(
                row,
                "upstreamEE",
            )

            downstream_es = coordinate(
                row,
                "downstreamES",
            )

            connections = [
                (
                    upstream_ee,
                    first_start,
                    "first_exon_upstream_junction",
                    "EXCISED_INTRON",
                ),
                (
                    first_end,
                    downstream_es,
                    "first_exon_downstream_junction",
                    "EXCISED_INTRON",
                ),
                (
                    upstream_ee,
                    second_start,
                    "second_exon_upstream_junction",
                    "EXCISED_INTRON",
                ),
                (
                    second_end,
                    downstream_es,
                    "second_exon_downstream_junction",
                    "EXCISED_INTRON",
                ),
            ]

        elif event_type in {
            "A3SS",
            "A5SS",
        }:
            long_start = coordinate(
                row,
                "longExonStart_0base",
            )

            long_end = coordinate(
                row,
                "longExonEnd",
            )

            short_start = coordinate(
                row,
                "shortES",
            )

            short_end = coordinate(
                row,
                "shortEE",
            )

            flanking_start = coordinate(
                row,
                "flankingES",
            )

            flanking_end = coordinate(
                row,
                "flankingEE",
            )

            variable_exon_is_upstream = (
                (
                    event_type == "A5SS"
                    and strand == "+"
                )
                or (
                    event_type == "A3SS"
                    and strand == "-"
                )
            )

            if variable_exon_is_upstream:
                connections = [
                    (
                        long_end,
                        flanking_start,
                        "long_isoform_junction",
                        "EXCISED_INTRON",
                    ),
                    (
                        short_end,
                        flanking_start,
                        "short_isoform_junction",
                        "EXCISED_INTRON",
                    ),
                ]

            else:
                connections = [
                    (
                        flanking_end,
                        long_start,
                        "long_isoform_junction",
                        "EXCISED_INTRON",
                    ),
                    (
                        flanking_end,
                        short_start,
                        "short_isoform_junction",
                        "EXCISED_INTRON",
                    ),
                ]

        else:
            upstream_ee = coordinate(
                row,
                "upstreamEE",
            )

            downstream_es = coordinate(
                row,
                "downstreamES",
            )

            connections = [
                (
                    upstream_ee,
                    downstream_es,
                    "retained_intron",
                    "RETAINED_INTRON",
                )
            ]

        for (
            start,
            end,
            role,
            feature_class,
        ) in connections:

            rows.append(
                {
                    "dataset_id": (
                        row[dataset_column]
                        if dataset_column is not None
                        else pd.NA
                    ),
                    "contrast": (
                        row[contrast_column]
                        if contrast_column is not None
                        else contrast
                    ),
                    "tool": "RMATS",
                    "tool_event_id": (
                        event_id
                    ),
                    "tool_group_id": (
                        event_id
                    ),
                    "connection_role": (
                        role
                    ),
                    "gene_id": (
                        row[gene_id_column]
                        if gene_id_column is not None
                        else pd.NA
                    ),
                    "gene_name": (
                        row[gene_name_column]
                        if gene_name_column is not None
                        else pd.NA
                    ),
                    "chromosome": (
                        row["chr"]
                    ),
                    "strand": (
                        strand
                    ),
                    "intron_start": min(
                        start,
                        end,
                    ),
                    "intron_end": max(
                        start,
                        end,
                    ),
                    "feature_class": (
                        feature_class
                    ),
                    "event_type": (
                        event_type
                    ),
                    "effect_value": (
                        effect
                    ),
                    "effect_metric": (
                        "dPSI_target_minus_reference"
                        if effect_column
                        == "dpsi_target_minus_reference"
                        else effect_column
                    ),
                    "event_significant": (
                        significant
                    ),
                    "connection_significant": (
                        significant
                    ),
                    "significance_value": (
                        fdr
                    ),
                    "significance_metric": (
                        "FDR"
                    ),
                    "source_file": (
                        row["source_file"]
                    ),
                    "source_row": (
                        row["source_row"]
                    ),
                }
            )

    return canonical_frame(
        rows
    )


def parse_rmats(
    paths: list[Path],
    contrast: str,
    *,
    fdr_threshold: float,
) -> pd.DataFrame:

    return pd.concat(
        [
            parse_rmats_file(
                path,
                contrast,
                fdr_threshold=fdr_threshold,
            )
            for path in paths
        ],
        ignore_index=True,
        sort=False,
    )


def parse_leafcutter(
    path: Path,
    contrast: str,
    *,
    fdr_threshold: float,
) -> pd.DataFrame:

    frame = read_tsv(
        path
    )

    cluster_column = find_column(
        frame,
        [
            "cluster",
            "cluster_id",
            "cluster_name",
        ],
        required=True,
        label="LeafCutter cluster",
        path=path,
    )

    intron_id_column = find_column(
        frame,
        [
            "intron_id",
            "intron",
            "junction_id",
        ],
        required=True,
        label="LeafCutter intron ID",
        path=path,
    )

    chromosome_column = find_column(
        frame,
        [
            "seqname",
            "chromosome",
            "chrom",
            "chr",
            "seqid",
            "seqname_gtf",
        ],
        required=True,
        label="LeafCutter chromosome",
        path=path,
    )

    start_column = find_column(
        frame,
        [
            "intron_start",
            "junction_start",
            "start",
        ],
        required=True,
        label="LeafCutter intron start",
        path=path,
    )

    end_column = find_column(
        frame,
        [
            "intron_end",
            "junction_end",
            "end",
        ],
        required=True,
        label="LeafCutter intron end",
        path=path,
    )

    strand_column = find_column(
        frame,
        [
            "strand",
            "gene_strand",
            "intron_strand",
        ],
        required=True,
        label="LeafCutter strand",
        path=path,
    )

    gene_id_column = find_column(
        frame,
        [
            "gene_id",
            "geneid",
        ],
    )

    gene_name_column = find_column(
        frame,
        [
            "gene_name",
            "gene_symbol",
            "genesymbol",
            "symbol",
            "gene",
        ],
    )

    effect_column = find_column(
        frame,
        [
            "deltapsi",
            "delta_psi",
            "dpsi",
            "delta_intron_usage",
            "effect_size",
        ],
    )

    adjusted_p_column = find_column(
        frame,
        [
            "cluster_p_adjust",
            "cluster_fdr",
            "fdr",
            "q_value",
            "qvalue",
            "padj",
        ],
    )

    cluster_flag_column = find_column(
        frame,
        [
            "cluster_significant",
            "significant",
            "is_significant",
        ],
    )

    intron_flag_column = find_column(
        frame,
        [
            "significant_intron",
            "intron_significant",
            "passes_effect_threshold",
        ],
    )

    contrast_column = find_column(
        frame,
        [
            "contrast",
            "comparison",
        ],
    )

    dataset_column = find_column(
        frame,
        [
            "dataset_id",
            "dataset",
        ],
    )

    starts = pd.to_numeric(
        frame[start_column],
        errors="coerce",
    )

    ends = pd.to_numeric(
        frame[end_column],
        errors="coerce",
    )

    adjusted_p = numeric_series(
        frame,
        adjusted_p_column,
    )

    event_significant = fill_significance(
        boolean_series(
            frame,
            cluster_flag_column,
        ),
        adjusted_p,
        mode="minimum",
        threshold=fdr_threshold,
    )

    connection_significant = boolean_series(
        frame,
        intron_flag_column,
    )

    connection_significant = (
        connection_significant.fillna(
            event_significant
        )
    )

    rows: list[dict[str, Any]] = []

    for index in frame.index:
        rows.append(
            {
                "dataset_id": (
                    frame.at[
                        index,
                        dataset_column,
                    ]
                    if dataset_column is not None
                    else pd.NA
                ),
                "contrast": (
                    frame.at[
                        index,
                        contrast_column,
                    ]
                    if contrast_column is not None
                    else contrast
                ),
                "tool": "LEAFCUTTER",
                "tool_event_id": (
                    frame.at[
                        index,
                        intron_id_column,
                    ]
                ),
                "tool_group_id": (
                    frame.at[
                        index,
                        cluster_column,
                    ]
                ),
                "connection_role": (
                    "cluster_intron"
                ),
                "gene_id": (
                    frame.at[
                        index,
                        gene_id_column,
                    ]
                    if gene_id_column is not None
                    else pd.NA
                ),
                "gene_name": (
                    frame.at[
                        index,
                        gene_name_column,
                    ]
                    if gene_name_column is not None
                    else pd.NA
                ),
                "chromosome": (
                    frame.at[
                        index,
                        chromosome_column,
                    ]
                ),
                "strand": (
                    frame.at[
                        index,
                        strand_column,
                    ]
                ),
                "intron_start": (
                    starts.at[index]
                ),
                "intron_end": (
                    ends.at[index]
                ),
                "feature_class": (
                    "EXCISED_INTRON"
                ),
                "event_type": (
                    "LEAFCUTTER_CLUSTER"
                ),
                "effect_value": (
                    frame.at[
                        index,
                        effect_column,
                    ]
                    if effect_column is not None
                    else pd.NA
                ),
                "effect_metric": (
                    effect_column
                    or pd.NA
                ),
                "event_significant": (
                    event_significant.at[
                        index
                    ]
                ),
                "connection_significant": (
                    connection_significant.at[
                        index
                    ]
                ),
                "significance_value": (
                    adjusted_p.at[index]
                ),
                "significance_metric": (
                    adjusted_p_column
                    or pd.NA
                ),
                "source_file": (
                    frame.at[
                        index,
                        "source_file",
                    ]
                ),
                "source_row": (
                    frame.at[
                        index,
                        "source_row",
                    ]
                ),
            }
        )

    return canonical_frame(
        rows
    )


def normalise_connections(
    frame: pd.DataFrame,
) -> pd.DataFrame:

    frame = frame.copy()

    for column in CANONICAL_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA

    text_columns = [
        "dataset_id",
        "contrast",
        "tool",
        "tool_event_id",
        "tool_group_id",
        "connection_role",
        "gene_name",
        "event_type",
        "effect_metric",
        "significance_metric",
        "source_file",
    ]

    for column in text_columns:
        frame[column] = clean_text(
            frame[column]
        )

    frame["gene_id"] = (
        clean_text(
            frame["gene_id"]
        )
        .str.replace(
            r"\.\d+$",
            "",
            regex=True,
        )
    )

    frame["tool"] = (
        frame["tool"]
        .str.upper()
    )

    chromosome = clean_text(
        frame["chromosome"]
    )

    frame["chromosome"] = (
        chromosome.map(
            lambda value: (
                value
                if (
                    pd.isna(value)
                    or str(value)
                    .lower()
                    .startswith("chr")
                )
                else f"chr{value}"
            )
        )
        .astype("string")
    )

    frame["strand"] = (
        clean_text(
            frame["strand"]
        )
        .replace(
            {
                "1": "+",
                "+1": "+",
                "-1": "-",
                "plus": "+",
                "minus": "-",
                "forward": "+",
                "reverse": "-",
            }
        )
    )

    frame["feature_class"] = (
        clean_text(
            frame["feature_class"]
        )
        .str.upper()
    )

    for column in [
        "intron_start",
        "intron_end",
        "source_row",
    ]:
        frame[column] = pd.to_numeric(
            frame[column],
            errors="coerce",
        ).astype("Int64")

    for column in [
        "effect_value",
        "significance_value",
    ]:
        frame[column] = pd.to_numeric(
            frame[column],
            errors="coerce",
        )

    for column in [
        "event_significant",
        "connection_significant",
    ]:
        frame[column] = (
            frame[column]
            .map(parse_boolean)
            .astype("boolean")
        )

    required_columns = [
        "tool",
        "tool_event_id",
        "tool_group_id",
        "chromosome",
        "strand",
        "intron_start",
        "intron_end",
        "feature_class",
    ]

    missing_required = (
        frame[required_columns]
        .isna()
        .any(axis=1)
    )

    if missing_required.any():
        raise ValueError(
            "Required canonical values are missing. "
            "Examples:\n"
            + frame.loc[
                missing_required,
                required_columns
                + [
                    "source_file",
                    "source_row",
                ],
            ]
            .head(10)
            .to_string(index=False)
        )

    invalid_strands = sorted(
        set(
            frame["strand"]
            .dropna()
        )
        - {
            "+",
            "-",
        }
    )

    if invalid_strands:
        raise ValueError(
            "Unsupported strand values: "
            + ", ".join(
                map(
                    str,
                    invalid_strands,
                )
            )
        )

    invalid_intervals = (
        frame["intron_start"]
        >= frame["intron_end"]
    )

    if invalid_intervals.any():
        raise ValueError(
            "Every connection must satisfy "
            "intron_start < intron_end. Examples:\n"
            + frame.loc[
                invalid_intervals,
                [
                    "tool",
                    "tool_event_id",
                    "chromosome",
                    "intron_start",
                    "intron_end",
                    "source_file",
                    "source_row",
                ],
            ]
            .head(10)
            .to_string(index=False)
        )

    frame["donor"] = (
        frame["intron_start"]
    )

    frame["acceptor"] = (
        frame["intron_end"]
    )

    negative = (
        frame["strand"]
        .eq("-")
    )

    frame.loc[
        negative,
        "donor",
    ] = frame.loc[
        negative,
        "intron_end",
    ]

    frame.loc[
        negative,
        "acceptor",
    ] = frame.loc[
        negative,
        "intron_start",
    ]

    frame["intron_locus_id"] = (
        frame["chromosome"].astype(str)
        + "|"
        + frame["strand"].astype(str)
        + "|"
        + frame["intron_start"].astype(str)
        + "|"
        + frame["intron_end"].astype(str)
    )

    frame["connection_id"] = (
        frame["intron_locus_id"]
        + "|"
        + frame["feature_class"].astype(str)
    )

    frame = (
        frame.drop_duplicates(
            subset=[
                "dataset_id",
                "contrast",
                "tool",
                "tool_event_id",
                "tool_group_id",
                "connection_role",
                "chromosome",
                "strand",
                "intron_start",
                "intron_end",
                "feature_class",
            ],
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    frame.insert(
        0,
        "evidence_id",
        [
            f"EVIDENCE_{index:08d}"
            for index in range(
                1,
                len(frame) + 1,
            )
        ],
    )

    frame["_scope"] = (
        frame["dataset_id"]
        .fillna("UNSPECIFIED")
        .astype(str)
        + "||"
        + frame["contrast"]
        .fillna("UNSPECIFIED")
        .astype(str)
    )

    frame["event_key"] = (
        frame["_scope"]
        + "||"
        + frame["tool"].astype(str)
        + "||"
        + frame["tool_group_id"].astype(str)
        + "||"
        + frame["chromosome"].astype(str)
        + "||"
        + frame["strand"].astype(str)
    )

    return frame


def join_unique(
    series: pd.Series,
    separator: str = "|",
) -> Any:

    values = sorted(
        {
            str(value)
            for value in series
            if (
                pd.notna(value)
                and str(value).strip()
            )
        }
    )

    return (
        separator.join(values)
        if values
        else pd.NA
    )


def summarise_boolean(
    series: pd.Series,
) -> Any:

    non_missing = (
        series.dropna()
    )

    if non_missing.empty:
        return pd.NA

    return bool(
        non_missing.astype(
            bool
        ).any()
    )


def is_true(
    value: Any,
) -> bool:

    return (
        False
        if pd.isna(value)
        else bool(value)
    )


def combine_delimited_values(
    values: Iterable[Any],
) -> Any:

    combined: set[str] = set()

    for value in values:
        if pd.isna(value):
            continue

        combined.update(
            part.strip()
            for part in str(value).split("|")
            if part.strip()
        )

    return (
        "|".join(
            sorted(combined)
        )
        if combined
        else pd.NA
    )


def build_event_regions(
    connections: pd.DataFrame,
) -> pd.DataFrame:

    group_columns = [
        "event_key",
        "_scope",
        "dataset_id",
        "contrast",
        "tool",
        "tool_group_id",
        "chromosome",
        "strand",
    ]

    rows: list[dict[str, Any]] = []

    for key, group in connections.groupby(
        group_columns,
        dropna=False,
        sort=False,
    ):
        row = dict(
            zip(
                group_columns,
                key,
            )
        )

        effects = pd.to_numeric(
            group["effect_value"],
            errors="coerce",
        ).dropna()

        if not effects.empty:
            maximum_effect_index = (
                effects.abs().idxmax()
            )

            max_abs_effect = float(
                effects.abs().max()
            )

            effect_at_max_abs = float(
                effects.loc[
                    maximum_effect_index
                ]
            )

        else:
            max_abs_effect = pd.NA
            effect_at_max_abs = pd.NA

        region_start = int(
            group["intron_start"].min()
        )

        region_end = int(
            group["intron_end"].max()
        )

        row.update(
            {
                "gene_ids": join_unique(
                    group["gene_id"]
                ),
                "gene_names": join_unique(
                    group["gene_name"]
                ),
                "region_start": (
                    region_start
                ),
                "region_end": (
                    region_end
                ),
                "region_length": (
                    region_end
                    - region_start
                ),
                "event_ids": join_unique(
                    group["tool_event_id"]
                ),
                "event_types": join_unique(
                    group["event_type"]
                ),
                "n_event_ids": (
                    group[
                        "tool_event_id"
                    ].nunique()
                ),
                "n_connections": (
                    group[
                        "connection_id"
                    ].nunique()
                ),
                "n_excised_introns": int(
                    group[
                        "feature_class"
                    ]
                    .eq(
                        "EXCISED_INTRON"
                    )
                    .sum()
                ),
                "n_retained_introns": int(
                    group[
                        "feature_class"
                    ]
                    .eq(
                        "RETAINED_INTRON"
                    )
                    .sum()
                ),
                "event_significant": (
                    summarise_boolean(
                        group[
                            "event_significant"
                        ]
                    )
                ),
                "n_significant_connections": int(
                    group[
                        "connection_significant"
                    ]
                    .fillna(False)
                    .sum()
                ),
                "effect_metrics": join_unique(
                    group[
                        "effect_metric"
                    ]
                ),
                "effect_min": (
                    effects.min()
                    if not effects.empty
                    else pd.NA
                ),
                "effect_max": (
                    effects.max()
                    if not effects.empty
                    else pd.NA
                ),
                "max_abs_effect": (
                    max_abs_effect
                ),
                "effect_at_max_abs": (
                    effect_at_max_abs
                ),
                "significance_metrics": join_unique(
                    group[
                        "significance_metric"
                    ]
                ),
                "source_files": join_unique(
                    group[
                        "source_file"
                    ]
                ),
            }
        )

        rows.append(
            row
        )

    regions = pd.DataFrame(
        rows
    )

    regions = regions.sort_values(
        [
            "dataset_id",
            "contrast",
            "chromosome",
            "strand",
            "region_start",
            "region_end",
            "tool",
            "tool_group_id",
        ],
        na_position="last",
        kind="mergesort",
    ).reset_index(
        drop=True
    )

    regions.insert(
        0,
        "tool_region_id",
        [
            f"TOOL_REGION_{index:08d}"
            for index in range(
                1,
                len(regions) + 1,
            )
        ],
    )

    return regions


def format_loci(
    values: set[tuple[Any, ...]],
) -> Any:

    if not values:
        return pd.NA

    return "|".join(
        ":".join(
            map(
                str,
                value,
            )
        )
        for value in sorted(values)
    )


def structural_support(
    connections_a: pd.DataFrame,
    connections_b: pd.DataFrame,
) -> dict[str, Any]:

    exact_a = {
        (
            row.chromosome,
            row.strand,
            int(row.intron_start),
            int(row.intron_end),
            row.feature_class,
        )
        for row in connections_a.itertuples(
            index=False
        )
    }

    exact_b = {
        (
            row.chromosome,
            row.strand,
            int(row.intron_start),
            int(row.intron_end),
            row.feature_class,
        )
        for row in connections_b.itertuples(
            index=False
        )
    }

    exact_connections = (
        exact_a
        & exact_b
    )

    excised_a = {
        (
            row.chromosome,
            row.strand,
            int(row.intron_start),
            int(row.intron_end),
        )
        for row in connections_a.loc[
            connections_a[
                "feature_class"
            ].eq(
                "EXCISED_INTRON"
            )
        ].itertuples(
            index=False
        )
    }

    excised_b = {
        (
            row.chromosome,
            row.strand,
            int(row.intron_start),
            int(row.intron_end),
        )
        for row in connections_b.loc[
            connections_b[
                "feature_class"
            ].eq(
                "EXCISED_INTRON"
            )
        ].itertuples(
            index=False
        )
    }

    retained_a = {
        (
            row.chromosome,
            row.strand,
            int(row.intron_start),
            int(row.intron_end),
        )
        for row in connections_a.loc[
            connections_a[
                "feature_class"
            ].eq(
                "RETAINED_INTRON"
            )
        ].itertuples(
            index=False
        )
    }

    retained_b = {
        (
            row.chromosome,
            row.strand,
            int(row.intron_start),
            int(row.intron_end),
        )
        for row in connections_b.loc[
            connections_b[
                "feature_class"
            ].eq(
                "RETAINED_INTRON"
            )
        ].itertuples(
            index=False
        )
    }

    same_intron_loci = (
        (
            excised_a
            & retained_b
        )
        |
        (
            retained_a
            & excised_b
        )
    )

    excised_frame_a = connections_a.loc[
        connections_a[
            "feature_class"
        ].eq(
            "EXCISED_INTRON"
        )
    ]

    excised_frame_b = connections_b.loc[
        connections_b[
            "feature_class"
        ].eq(
            "EXCISED_INTRON"
        )
    ]

    donor_a = set(
        pd.to_numeric(
            excised_frame_a[
                "donor"
            ],
            errors="coerce",
        )
        .dropna()
        .astype(int)
    )

    donor_b = set(
        pd.to_numeric(
            excised_frame_b[
                "donor"
            ],
            errors="coerce",
        )
        .dropna()
        .astype(int)
    )

    acceptor_a = set(
        pd.to_numeric(
            excised_frame_a[
                "acceptor"
            ],
            errors="coerce",
        )
        .dropna()
        .astype(int)
    )

    acceptor_b = set(
        pd.to_numeric(
            excised_frame_b[
                "acceptor"
            ],
            errors="coerce",
        )
        .dropna()
        .astype(int)
    )

    shared_donors = (
        donor_a
        & donor_b
    )

    shared_acceptors = (
        acceptor_a
        & acceptor_b
    )

    if exact_connections:
        support_class = (
            "EXACT_CONNECTION"
        )

    elif same_intron_loci:
        support_class = (
            "SAME_INTRON_LOCUS"
        )

    elif (
        shared_donors
        or shared_acceptors
    ):
        support_class = (
            "SHARED_SPLICE_SITE"
        )

    else:
        support_class = (
            "REGIONAL_OVERLAP"
        )

    return {
        "support_class": (
            support_class
        ),
        "support_priority": (
            SUPPORT_PRIORITY[
                support_class
            ]
        ),
        "n_exact_connections": (
            len(
                exact_connections
            )
        ),
        "n_same_intron_loci": (
            len(
                same_intron_loci
            )
        ),
        "n_shared_donor_sites": (
            len(
                shared_donors
            )
        ),
        "n_shared_acceptor_sites": (
            len(
                shared_acceptors
            )
        ),
        "exact_connections": (
            format_loci(
                exact_connections
            )
        ),
        "same_intron_loci": (
            format_loci(
                same_intron_loci
            )
        ),
        "shared_donor_sites": (
            "|".join(
                map(
                    str,
                    sorted(
                        shared_donors
                    ),
                )
            )
            if shared_donors
            else pd.NA
        ),
        "shared_acceptor_sites": (
            "|".join(
                map(
                    str,
                    sorted(
                        shared_acceptors
                    ),
                )
            )
            if shared_acceptors
            else pd.NA
        ),
    }


def compare_event_regions(
    regions: pd.DataFrame,
    connections: pd.DataFrame,
    *,
    overlap_threshold: float,
) -> pd.DataFrame:

    connections_by_event = {
        event_key: group.copy()
        for event_key, group
        in connections.groupby(
            "event_key",
            sort=False,
        )
    }

    rows: list[dict[str, Any]] = []

    partition_columns = [
        "_scope",
        "chromosome",
        "strand",
    ]

    side_columns = [
        "tool_region_id",
        "tool",
        "tool_group_id",
        "gene_ids",
        "gene_names",
        "region_start",
        "region_end",
        "region_length",
        "event_ids",
        "event_types",
        "n_event_ids",
        "n_connections",
        "event_significant",
        "n_significant_connections",
        "effect_metrics",
        "effect_min",
        "effect_max",
        "max_abs_effect",
        "effect_at_max_abs",
        "significance_metrics",
    ]

    for _, partition in regions.groupby(
        partition_columns,
        sort=False,
        dropna=False,
    ):
        ordered = partition.sort_values(
            [
                "region_start",
                "region_end",
            ],
            kind="mergesort",
        )

        active: list[int] = []

        for current_index in ordered.index:
            current_start = int(
                regions.at[
                    current_index,
                    "region_start",
                ]
            )

            active = [
                index
                for index in active
                if int(
                    regions.at[
                        index,
                        "region_end",
                    ]
                )
                > current_start
            ]

            for other_index in active:
                left = regions.loc[
                    other_index
                ]

                right = regions.loc[
                    current_index
                ]

                if (
                    left["tool"]
                    == right["tool"]
                ):
                    continue

                if (
                    str(left["tool"]),
                    str(
                        left[
                            "tool_group_id"
                        ]
                    ),
                ) > (
                    str(right["tool"]),
                    str(
                        right[
                            "tool_group_id"
                        ]
                    ),
                ):
                    left, right = (
                        right,
                        left,
                    )

                overlap_start = max(
                    int(
                        left[
                            "region_start"
                        ]
                    ),
                    int(
                        right[
                            "region_start"
                        ]
                    ),
                )

                overlap_end = min(
                    int(
                        left[
                            "region_end"
                        ]
                    ),
                    int(
                        right[
                            "region_end"
                        ]
                    ),
                )

                overlap_bp = max(
                    0,
                    overlap_end
                    - overlap_start,
                )

                if overlap_bp == 0:
                    continue

                fraction_a = (
                    overlap_bp
                    / int(
                        left[
                            "region_length"
                        ]
                    )
                )

                fraction_b = (
                    overlap_bp
                    / int(
                        right[
                            "region_length"
                        ]
                    )
                )

                reciprocal_overlap = min(
                    fraction_a,
                    fraction_b,
                )

                passes_threshold = (
                    reciprocal_overlap
                    >= overlap_threshold
                )

                significant_a = is_true(
                    left[
                        "event_significant"
                    ]
                )

                significant_b = is_true(
                    right[
                        "event_significant"
                    ]
                )

                both_significant = (
                    significant_a
                    and significant_b
                )

                support = structural_support(
                    connections_by_event[
                        left[
                            "event_key"
                        ]
                    ],
                    connections_by_event[
                        right[
                            "event_key"
                        ]
                    ],
                )

                gene_ids = (
                    combine_delimited_values(
                        [
                            left[
                                "gene_ids"
                            ],
                            right[
                                "gene_ids"
                            ],
                        ]
                    )
                )

                gene_names = (
                    combine_delimited_values(
                        [
                            left[
                                "gene_names"
                            ],
                            right[
                                "gene_names"
                            ],
                        ]
                    )
                )

                row: dict[str, Any] = {
                    "dataset_id": (
                        left[
                            "dataset_id"
                        ]
                        if pd.notna(
                            left[
                                "dataset_id"
                            ]
                        )
                        else right[
                            "dataset_id"
                        ]
                    ),
                    "contrast": (
                        left[
                            "contrast"
                        ]
                        if pd.notna(
                            left[
                                "contrast"
                            ]
                        )
                        else right[
                            "contrast"
                        ]
                    ),
                    "chromosome": (
                        left[
                            "chromosome"
                        ]
                    ),
                    "strand": (
                        left[
                            "strand"
                        ]
                    ),
                    "gene_ids": (
                        gene_ids
                    ),
                    "gene_names": (
                        gene_names
                    ),
                    "overlap_start": (
                        overlap_start
                    ),
                    "overlap_end": (
                        overlap_end
                    ),
                    "overlap_bp": (
                        overlap_bp
                    ),
                    "overlap_fraction_a": (
                        fraction_a
                    ),
                    "overlap_fraction_b": (
                        fraction_b
                    ),
                    "reciprocal_overlap": (
                        reciprocal_overlap
                    ),
                    "overlap_threshold": (
                        overlap_threshold
                    ),
                    "passes_overlap_threshold": (
                        passes_threshold
                    ),
                    "both_events_significant": (
                        both_significant
                    ),
                    "significant_overlap_hit": (
                        passes_threshold
                        and both_significant
                    ),
                    **support,
                    "event_key_a": (
                        left[
                            "event_key"
                        ]
                    ),
                    "event_key_b": (
                        right[
                            "event_key"
                        ]
                    ),
                }

                for suffix, region in [
                    (
                        "a",
                        left,
                    ),
                    (
                        "b",
                        right,
                    ),
                ]:
                    for column in side_columns:
                        row[
                            f"{column}_{suffix}"
                        ] = region[
                            column
                        ]

                rows.append(
                    row
                )

            active.append(
                current_index
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                "dataset_id",
                "contrast",
                "chromosome",
                "strand",
                "gene_ids",
                "gene_names",
                "overlap_start",
                "overlap_end",
                "overlap_bp",
                "overlap_fraction_a",
                "overlap_fraction_b",
                "reciprocal_overlap",
                "passes_overlap_threshold",
                "both_events_significant",
                "significant_overlap_hit",
                "support_class",
            ]
        )

    return (
        pd.DataFrame(
            rows
        )
        .sort_values(
            [
                "significant_overlap_hit",
                "reciprocal_overlap",
                "support_priority",
                "chromosome",
                "overlap_start",
                "tool_a",
                "tool_b",
            ],
            ascending=[
                False,
                False,
                True,
                True,
                True,
                True,
                True,
            ],
            na_position="last",
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )


def build_multi_tool_summary(
    significant_overlaps: pd.DataFrame,
    regions: pd.DataFrame,
) -> pd.DataFrame:

    empty_columns = [
        "multi_tool_region_id",
        "dataset_id",
        "contrast",
        "chromosome",
        "strand",
        "region_start",
        "region_end",
        "region_length",
        "gene_ids",
        "gene_names",
        "tools",
        "n_tools",
        "significant_tools",
        "n_tools_significant",
        "n_tool_event_regions",
        "event_groups",
        "event_types",
        "member_event_keys",
        "strongest_support",
        "support_classes",
        "n_significant_overlap_pairs",
        "minimum_reciprocal_overlap",
        "maximum_reciprocal_overlap",
        "n_exact_connections",
        "n_same_intron_loci",
        "n_shared_donor_sites",
        "n_shared_acceptor_sites",
        "effect_metrics",
        "max_abs_effect",
        "effect_at_max_abs",
        "significance_metrics",
        "source_files",
    ]

    if significant_overlaps.empty:
        return pd.DataFrame(
            columns=empty_columns
        )

    nodes = (
        set(
            significant_overlaps[
                "event_key_a"
            ]
        )
        |
        set(
            significant_overlaps[
                "event_key_b"
            ]
        )
    )

    union_find = UnionFind(
        nodes
    )

    for row in (
        significant_overlaps
        .itertuples(
            index=False
        )
    ):
        union_find.union(
            row.event_key_a,
            row.event_key_b,
        )

    component_members: dict[
        str,
        set[str],
    ] = defaultdict(
        set
    )

    for node in nodes:
        component_members[
            union_find.find(node)
        ].add(
            node
        )

    region_lookup = regions.set_index(
        "event_key",
        drop=False,
    )

    rows: list[dict[str, Any]] = []

    for component, members in (
        component_members.items()
    ):
        member_regions = region_lookup.loc[
            sorted(members)
        ]

        component_edges = (
            significant_overlaps.loc[
                significant_overlaps[
                    "event_key_a"
                ].isin(members)
                &
                significant_overlaps[
                    "event_key_b"
                ].isin(members)
            ]
        )

        support_classes = sorted(
            set(
                component_edges[
                    "support_class"
                ]
            ),
            key=lambda value: (
                SUPPORT_PRIORITY[
                    value
                ]
            ),
        )

        strongest_support = (
            support_classes[0]
        )

        significant_member_mask = (
            member_regions[
                "event_significant"
            ]
            .fillna(False)
            .astype(bool)
        )

        significant_tools = sorted(
            set(
                member_regions.loc[
                    significant_member_mask,
                    "tool",
                ]
                .dropna()
                .astype(str)
            )
        )

        event_groups = "|".join(
            sorted(
                f"{row.tool}:"
                f"{row.tool_group_id}"
                for row in (
                    member_regions
                    .itertuples(
                        index=False
                    )
                )
            )
        )

        region_start = int(
            member_regions[
                "region_start"
            ].min()
        )

        region_end = int(
            member_regions[
                "region_end"
            ].max()
        )

        member_effects = pd.to_numeric(
            member_regions[
                "effect_at_max_abs"
            ],
            errors="coerce",
        ).dropna()

        if not member_effects.empty:
            maximum_effect_index = (
                member_effects
                .abs()
                .idxmax()
            )

            max_abs_effect = float(
                member_effects
                .abs()
                .max()
            )

            effect_at_max_abs = float(
                member_effects.loc[
                    maximum_effect_index
                ]
            )

        else:
            max_abs_effect = pd.NA
            effect_at_max_abs = pd.NA

        rows.append(
            {
                "component_key": (
                    component
                ),
                "dataset_id": join_unique(
                    member_regions[
                        "dataset_id"
                    ]
                ),
                "contrast": join_unique(
                    member_regions[
                        "contrast"
                    ]
                ),
                "chromosome": join_unique(
                    member_regions[
                        "chromosome"
                    ]
                ),
                "strand": join_unique(
                    member_regions[
                        "strand"
                    ]
                ),
                "region_start": (
                    region_start
                ),
                "region_end": (
                    region_end
                ),
                "region_length": (
                    region_end
                    - region_start
                ),
                "gene_ids": (
                    combine_delimited_values(
                        member_regions[
                            "gene_ids"
                        ]
                    )
                ),
                "gene_names": (
                    combine_delimited_values(
                        member_regions[
                            "gene_names"
                        ]
                    )
                ),
                "tools": (
                    combine_delimited_values(
                        member_regions[
                            "tool"
                        ]
                    )
                ),
                "n_tools": (
                    member_regions[
                        "tool"
                    ].nunique()
                ),
                "significant_tools": (
                    "|".join(
                        significant_tools
                    )
                    if significant_tools
                    else pd.NA
                ),
                "n_tools_significant": (
                    len(
                        significant_tools
                    )
                ),
                "n_tool_event_regions": (
                    len(
                        member_regions
                    )
                ),
                "event_groups": (
                    event_groups
                ),
                "event_types": (
                    combine_delimited_values(
                        member_regions[
                            "event_types"
                        ]
                    )
                ),
                "member_event_keys": (
                    ";".join(
                        sorted(
                            members
                        )
                    )
                ),
                "strongest_support": (
                    strongest_support
                ),
                "support_priority": (
                    SUPPORT_PRIORITY[
                        strongest_support
                    ]
                ),
                "support_classes": (
                    "|".join(
                        support_classes
                    )
                ),
                "n_significant_overlap_pairs": (
                    len(
                        component_edges
                    )
                ),
                "minimum_reciprocal_overlap": (
                    component_edges[
                        "reciprocal_overlap"
                    ].min()
                ),
                "maximum_reciprocal_overlap": (
                    component_edges[
                        "reciprocal_overlap"
                    ].max()
                ),
                "n_exact_connections": int(
                    component_edges[
                        "n_exact_connections"
                    ].sum()
                ),
                "n_same_intron_loci": int(
                    component_edges[
                        "n_same_intron_loci"
                    ].sum()
                ),
                "n_shared_donor_sites": int(
                    component_edges[
                        "n_shared_donor_sites"
                    ].sum()
                ),
                "n_shared_acceptor_sites": int(
                    component_edges[
                        "n_shared_acceptor_sites"
                    ].sum()
                ),
                "effect_metrics": (
                    combine_delimited_values(
                        member_regions[
                            "effect_metrics"
                        ]
                    )
                ),
                "max_abs_effect": (
                    max_abs_effect
                ),
                "effect_at_max_abs": (
                    effect_at_max_abs
                ),
                "significance_metrics": (
                    combine_delimited_values(
                        member_regions[
                            "significance_metrics"
                        ]
                    )
                ),
                "source_files": (
                    combine_delimited_values(
                        member_regions[
                            "source_files"
                        ]
                    )
                ),
            }
        )

    summary = (
        pd.DataFrame(
            rows
        )
        .sort_values(
            [
                "n_tools_significant",
                "n_tools",
                "support_priority",
                "minimum_reciprocal_overlap",
                "max_abs_effect",
                "chromosome",
                "region_start",
            ],
            ascending=[
                False,
                False,
                True,
                False,
                False,
                True,
                True,
            ],
            na_position="last",
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )

    summary.insert(
        0,
        "multi_tool_region_id",
        [
            f"MULTI_TOOL_REGION_{index:06d}"
            for index in range(
                1,
                len(summary) + 1,
            )
        ],
    )

    return summary.drop(
        columns=[
            "component_key",
            "support_priority",
        ]
    )


def calculate_reciprocal_overlap(
    start_a: int,
    end_a: int,
    start_b: int,
    end_b: int,
) -> float:

    overlap_start = max(
        int(start_a),
        int(start_b),
    )

    overlap_end = min(
        int(end_a),
        int(end_b),
    )

    overlap_bp = max(
        0,
        overlap_end
        - overlap_start,
    )

    if overlap_bp == 0:
        return 0.0

    length_a = (
        int(end_a)
        - int(start_a)
    )

    length_b = (
        int(end_b)
        - int(start_b)
    )

    if (
        length_a <= 0
        or length_b <= 0
    ):
        return 0.0

    return min(
        overlap_bp / length_a,
        overlap_bp / length_b,
    )


def collapse_single_tool_regions(
    regions: pd.DataFrame,
    *,
    overlap_threshold: float,
) -> pd.DataFrame:
    """
    Collapse redundant significant regions from the same tool.

    A region joins a group only when it reciprocally overlaps that group's
    representative region by at least overlap_threshold. This avoids
    transitive chaining across a long locus.
    """

    if regions.empty:
        return pd.DataFrame()

    work = regions.copy()

    work["_effect_sort"] = (
        pd.to_numeric(
            work[
                "max_abs_effect"
            ],
            errors="coerce",
        )
        .fillna(-1.0)
    )

    rows: list[dict[str, Any]] = []

    partition_columns = [
        "_scope",
        "chromosome",
        "strand",
        "tool",
    ]

    for _, partition in work.groupby(
        partition_columns,
        sort=False,
        dropna=False,
    ):
        remaining = list(
            partition.sort_values(
                [
                    "_effect_sort",
                    "region_length",
                    "region_start",
                    "region_end",
                ],
                ascending=[
                    False,
                    True,
                    True,
                    True,
                ],
                kind="mergesort",
            ).index
        )

        while remaining:
            representative_index = (
                remaining[0]
            )

            representative = work.loc[
                representative_index
            ]

            member_indices: list[int] = []
            unassigned_indices: list[int] = []

            for candidate_index in remaining:
                candidate = work.loc[
                    candidate_index
                ]

                reciprocal_overlap = (
                    calculate_reciprocal_overlap(
                        representative[
                            "region_start"
                        ],
                        representative[
                            "region_end"
                        ],
                        candidate[
                            "region_start"
                        ],
                        candidate[
                            "region_end"
                        ],
                    )
                )

                if (
                    reciprocal_overlap
                    >= overlap_threshold
                ):
                    member_indices.append(
                        candidate_index
                    )

                else:
                    unassigned_indices.append(
                        candidate_index
                    )

            remaining = (
                unassigned_indices
            )

            members = work.loc[
                member_indices
            ]

            region_start = int(
                members[
                    "region_start"
                ].min()
            )

            region_end = int(
                members[
                    "region_end"
                ].max()
            )

            member_effects = pd.to_numeric(
                members[
                    "effect_at_max_abs"
                ],
                errors="coerce",
            ).dropna()

            if not member_effects.empty:
                maximum_effect_index = (
                    member_effects
                    .abs()
                    .idxmax()
                )

                max_abs_effect = float(
                    member_effects
                    .abs()
                    .max()
                )

                effect_at_max_abs = float(
                    member_effects.loc[
                        maximum_effect_index
                    ]
                )

            else:
                max_abs_effect = pd.NA
                effect_at_max_abs = pd.NA

            tool = str(
                representative[
                    "tool"
                ]
            )

            rows.append(
                {
                    "ranking_class": (
                        "SINGLE_TOOL_SIGNIFICANT"
                    ),
                    "dataset_id": join_unique(
                        members[
                            "dataset_id"
                        ]
                    ),
                    "contrast": join_unique(
                        members[
                            "contrast"
                        ]
                    ),
                    "chromosome": (
                        representative[
                            "chromosome"
                        ]
                    ),
                    "strand": (
                        representative[
                            "strand"
                        ]
                    ),
                    "region_start": (
                        region_start
                    ),
                    "region_end": (
                        region_end
                    ),
                    "region_length": (
                        region_end
                        - region_start
                    ),
                    "gene_ids": (
                        combine_delimited_values(
                            members[
                                "gene_ids"
                            ]
                        )
                    ),
                    "gene_names": (
                        combine_delimited_values(
                            members[
                                "gene_names"
                            ]
                        )
                    ),
                    "tools_detected": (
                        tool
                    ),
                    "n_tools_detected": 1,
                    "significant_tools": (
                        tool
                    ),
                    "n_tools_significant": 1,
                    "strongest_support": (
                        "SINGLE_TOOL"
                    ),
                    "support_classes": (
                        "SINGLE_TOOL"
                    ),
                    "minimum_reciprocal_overlap": (
                        pd.NA
                    ),
                    "maximum_reciprocal_overlap": (
                        pd.NA
                    ),
                    "n_exact_connections": 0,
                    "n_same_intron_loci": 0,
                    "n_shared_donor_sites": 0,
                    "n_shared_acceptor_sites": 0,
                    "n_significant_overlap_pairs": 0,
                    "n_tool_event_regions": (
                        len(
                            members
                        )
                    ),
                    "event_groups": (
                        join_unique(
                            members[
                                "tool_group_id"
                            ]
                        )
                    ),
                    "event_types": (
                        combine_delimited_values(
                            members[
                                "event_types"
                            ]
                        )
                    ),
                    "member_event_keys": (
                        ";".join(
                            sorted(
                                members[
                                    "event_key"
                                ]
                                .dropna()
                                .astype(str)
                            )
                        )
                    ),
                    "effect_metrics": (
                        combine_delimited_values(
                            members[
                                "effect_metrics"
                            ]
                        )
                    ),
                    "max_abs_effect": (
                        max_abs_effect
                    ),
                    "effect_at_max_abs": (
                        effect_at_max_abs
                    ),
                    "significance_metrics": (
                        combine_delimited_values(
                            members[
                                "significance_metrics"
                            ]
                        )
                    ),
                    "source_files": (
                        combine_delimited_values(
                            members[
                                "source_files"
                            ]
                        )
                    ),
                }
            )

    return pd.DataFrame(
        rows
    )


def build_ranked_splicing_events(
    multi_tool_summary: pd.DataFrame,
    significant_overlaps: pd.DataFrame,
    event_regions: pd.DataFrame,
    *,
    overlap_threshold: float,
    sashimi_top_n: int,
) -> pd.DataFrame:
    """
    Rank significant regions for reporting and sashimi plotting.

    Ranking order:

      1. Significant multi-tool regions.
      2. Remaining significant MAJIQ regions.
      3. Remaining significant LeafCutter regions.
      4. Remaining significant rMATS regions.

    Within each single-tool category, larger absolute effects rank first.
    Native statistical measures are retained but are not combined into a
    cross-tool numeric score.
    """

    multi_tool_rows: list[
        dict[str, Any]
    ] = []

    for row in (
        multi_tool_summary
        .to_dict(
            orient="records"
        )
    ):
        multi_tool_rows.append(
            {
                "ranking_class": (
                    "MULTI_TOOL_SIGNIFICANT"
                ),
                "dataset_id": (
                    row.get(
                        "dataset_id",
                        pd.NA,
                    )
                ),
                "contrast": (
                    row.get(
                        "contrast",
                        pd.NA,
                    )
                ),
                "chromosome": (
                    row[
                        "chromosome"
                    ]
                ),
                "strand": (
                    row[
                        "strand"
                    ]
                ),
                "region_start": (
                    row[
                        "region_start"
                    ]
                ),
                "region_end": (
                    row[
                        "region_end"
                    ]
                ),
                "region_length": (
                    row[
                        "region_length"
                    ]
                ),
                "gene_ids": (
                    row.get(
                        "gene_ids",
                        pd.NA,
                    )
                ),
                "gene_names": (
                    row.get(
                        "gene_names",
                        pd.NA,
                    )
                ),
                "tools_detected": (
                    row[
                        "tools"
                    ]
                ),
                "n_tools_detected": (
                    row[
                        "n_tools"
                    ]
                ),
                "significant_tools": (
                    row.get(
                        "significant_tools",
                        row[
                            "tools"
                        ],
                    )
                ),
                "n_tools_significant": (
                    row.get(
                        "n_tools_significant",
                        row[
                            "n_tools"
                        ],
                    )
                ),
                "strongest_support": (
                    row[
                        "strongest_support"
                    ]
                ),
                "support_classes": (
                    row[
                        "support_classes"
                    ]
                ),
                "minimum_reciprocal_overlap": (
                    row[
                        "minimum_reciprocal_overlap"
                    ]
                ),
                "maximum_reciprocal_overlap": (
                    row[
                        "maximum_reciprocal_overlap"
                    ]
                ),
                "n_exact_connections": (
                    row[
                        "n_exact_connections"
                    ]
                ),
                "n_same_intron_loci": (
                    row[
                        "n_same_intron_loci"
                    ]
                ),
                "n_shared_donor_sites": (
                    row[
                        "n_shared_donor_sites"
                    ]
                ),
                "n_shared_acceptor_sites": (
                    row[
                        "n_shared_acceptor_sites"
                    ]
                ),
                "n_significant_overlap_pairs": (
                    row[
                        "n_significant_overlap_pairs"
                    ]
                ),
                "n_tool_event_regions": (
                    row[
                        "n_tool_event_regions"
                    ]
                ),
                "event_groups": (
                    row[
                        "event_groups"
                    ]
                ),
                "event_types": (
                    row[
                        "event_types"
                    ]
                ),
                "member_event_keys": (
                    row[
                        "member_event_keys"
                    ]
                ),
                "effect_metrics": (
                    row.get(
                        "effect_metrics",
                        pd.NA,
                    )
                ),
                "max_abs_effect": (
                    row.get(
                        "max_abs_effect",
                        pd.NA,
                    )
                ),
                "effect_at_max_abs": (
                    row.get(
                        "effect_at_max_abs",
                        pd.NA,
                    )
                ),
                "significance_metrics": (
                    row.get(
                        "significance_metrics",
                        pd.NA,
                    )
                ),
                "source_files": (
                    row.get(
                        "source_files",
                        pd.NA,
                    )
                ),
            }
        )

    multi_tool_ranked = pd.DataFrame(
        multi_tool_rows
    )

    covered_event_keys: set[str] = set()

    if not significant_overlaps.empty:
        covered_event_keys.update(
            significant_overlaps[
                "event_key_a"
            ]
            .dropna()
            .astype(str)
        )

        covered_event_keys.update(
            significant_overlaps[
                "event_key_b"
            ]
            .dropna()
            .astype(str)
        )

    significant_mask = (
        event_regions[
            "event_significant"
        ]
        .fillna(False)
        .astype(bool)
    )

    uncovered_mask = (
        ~event_regions[
            "event_key"
        ]
        .astype(str)
        .isin(
            covered_event_keys
        )
    )

    single_tool_regions = (
        event_regions.loc[
            significant_mask
            & uncovered_mask
        ]
        .copy()
    )

    single_tool_ranked = (
        collapse_single_tool_regions(
            single_tool_regions,
            overlap_threshold=(
                overlap_threshold
            ),
        )
    )

    ranked = pd.concat(
        [
            multi_tool_ranked,
            single_tool_ranked,
        ],
        ignore_index=True,
        sort=False,
    )

    ordered_columns = [
        "rank",
        "default_sashimi_plot",
        "ranking_class",
        "gene_ids",
        "gene_names",
        "chromosome",
        "strand",
        "region_start",
        "region_end",
        "region_length",
        "sashimi_region",
        "tools_detected",
        "n_tools_detected",
        "significant_tools",
        "n_tools_significant",
        "strongest_support",
        "support_classes",
        "minimum_reciprocal_overlap",
        "maximum_reciprocal_overlap",
        "n_exact_connections",
        "n_same_intron_loci",
        "n_shared_donor_sites",
        "n_shared_acceptor_sites",
        "n_significant_overlap_pairs",
        "n_tool_event_regions",
        "event_groups",
        "event_types",
        "member_event_keys",
        "max_abs_effect",
        "effect_at_max_abs",
        "effect_metrics",
        "significance_metrics",
        "source_files",
        "dataset_id",
        "contrast",
    ]

    if ranked.empty:
        return pd.DataFrame(
            columns=ordered_columns
        )

    ranking_class_priority = {
        "MULTI_TOOL_SIGNIFICANT": 1,
        "SINGLE_TOOL_SIGNIFICANT": 2,
    }

    single_tool_priority = {
        "MAJIQ": 1,
        "LEAFCUTTER": 2,
        "RMATS": 3,
    }

    support_priority = {
        "EXACT_CONNECTION": 1,
        "SAME_INTRON_LOCUS": 2,
        "SHARED_SPLICE_SITE": 3,
        "REGIONAL_OVERLAP": 4,
        "SINGLE_TOOL": 5,
    }

    ranked[
        "_ranking_class_priority"
    ] = (
        ranked[
            "ranking_class"
        ]
        .map(
            ranking_class_priority
        )
        .fillna(99)
    )

    def determine_single_tool_priority(
        row: pd.Series,
    ) -> int:
        if (
            row[
                "ranking_class"
            ]
            != "SINGLE_TOOL_SIGNIFICANT"
        ):
            return 0

        tool = str(
            row.get(
                "significant_tools",
                row.get(
                    "tools_detected",
                    "",
                ),
            )
        ).strip().upper()

        return single_tool_priority.get(
            tool,
            99,
        )

    ranked[
        "_single_tool_priority"
    ] = ranked.apply(
        determine_single_tool_priority,
        axis=1,
    )

    ranked[
        "_support_priority"
    ] = (
        ranked[
            "strongest_support"
        ]
        .map(
            support_priority
        )
        .fillna(99)
    )

    ranked[
        "_overlap_sort"
    ] = (
        pd.to_numeric(
            ranked[
                "minimum_reciprocal_overlap"
            ],
            errors="coerce",
        )
        .fillna(-1.0)
    )

    ranked[
        "_effect_sort"
    ] = (
        pd.to_numeric(
            ranked[
                "max_abs_effect"
            ],
            errors="coerce",
        )
        .fillna(-1.0)
    )

    ranked[
        "n_tools_significant"
    ] = (
        pd.to_numeric(
            ranked[
                "n_tools_significant"
            ],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )

    ranked[
        "n_tools_detected"
    ] = (
        pd.to_numeric(
            ranked[
                "n_tools_detected"
            ],
            errors="coerce",
        )
        .fillna(0)
        .astype(int)
    )

    ranked = (
        ranked.sort_values(
            [
                "_ranking_class_priority",
                "_single_tool_priority",
                "n_tools_significant",
                "n_tools_detected",
                "_support_priority",
                "_overlap_sort",
                "_effect_sort",
                "gene_names",
                "chromosome",
                "region_start",
                "region_end",
            ],
            ascending=[
                True,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                True,
                True,
                True,
            ],
            na_position="last",
            kind="mergesort",
        )
        .reset_index(
            drop=True
        )
    )

    ranked.insert(
        0,
        "rank",
        range(
            1,
            len(ranked) + 1,
        ),
    )

    ranked[
        "default_sashimi_plot"
    ] = (
        ranked[
            "rank"
        ]
        <= int(
            sashimi_top_n
        )
    )

    ranked[
        "sashimi_region"
    ] = (
        ranked[
            "chromosome"
        ].astype(str)
        + ":"
        + ranked[
            "region_start"
        ].astype(
            "Int64"
        ).astype(str)
        + "-"
        + ranked[
            "region_end"
        ].astype(
            "Int64"
        ).astype(str)
    )

    ranked = ranked.drop(
        columns=[
            "_ranking_class_priority",
            "_single_tool_priority",
            "_support_priority",
            "_overlap_sort",
            "_effect_sort",
        ],
        errors="ignore",
    )

    for column in ordered_columns:
        if column not in ranked.columns:
            ranked[
                column
            ] = pd.NA

    return ranked[
        ordered_columns
    ]


def write_excel_sheet(
    writer: pd.ExcelWriter,
    frame: pd.DataFrame,
    sheet_name: str,
) -> None:

    if len(frame) <= EXCEL_MAX_DATA_ROWS:
        frame.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )

    else:
        pd.DataFrame(
            {
                "message": [
                    (
                        f"{sheet_name} contains "
                        f"{len(frame):,} rows and exceeds "
                        f"the Excel worksheet limit. "
                        f"Use the TSV output."
                    )
                ]
            }
        ).to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )


def write_workbook(
    path: Path,
    ranked_events: pd.DataFrame,
    multi_tool_summary: pd.DataFrame,
    significant_overlaps: pd.DataFrame,
    all_overlaps: pd.DataFrame,
    event_regions: pd.DataFrame,
    connections: pd.DataFrame,
) -> None:

    with pd.ExcelWriter(
        path,
        engine="openpyxl",
    ) as writer:

        write_excel_sheet(
            writer,
            ranked_events,
            "Ranked splicing events",
        )

        write_excel_sheet(
            writer,
            multi_tool_summary,
            "Multi-tool regions",
        )

        write_excel_sheet(
            writer,
            significant_overlaps,
            "Significant overlaps",
        )

        write_excel_sheet(
            writer,
            all_overlaps,
            "All overlaps",
        )

        write_excel_sheet(
            writer,
            event_regions.drop(
                columns=[
                    "_scope",
                    "event_key",
                ],
                errors="ignore",
            ),
            "Tool event regions",
        )

        write_excel_sheet(
            writer,
            connections.drop(
                columns=[
                    "_scope",
                    "event_key",
                ],
                errors="ignore",
            ),
            "All connections",
        )

        for worksheet in (
            writer.book.worksheets
        ):
            worksheet.freeze_panes = (
                "A2"
            )

            worksheet.auto_filter.ref = (
                worksheet.dimensions
            )


def main() -> int:
    args = parse_arguments()

    connections = normalise_connections(
        pd.concat(
            [
                parse_majiq(
                    args.majiq_junctions,
                    args.prefix,
                    probability_threshold=(
                        args.majiq_probability_threshold
                    ),
                ),
                parse_rmats(
                    args.rmats_jc,
                    args.prefix,
                    fdr_threshold=(
                        args.fdr_threshold
                    ),
                ),
                parse_leafcutter(
                    args.leafcutter_introns,
                    args.prefix,
                    fdr_threshold=(
                        args.fdr_threshold
                    ),
                ),
            ],
            ignore_index=True,
            sort=False,
        )
    )

    event_regions = build_event_regions(
        connections
    )

    all_overlaps = compare_event_regions(
        event_regions,
        connections,
        overlap_threshold=(
            args.overlap_threshold
        ),
    )

    significant_overlaps = (
        all_overlaps.loc[
            all_overlaps[
                "significant_overlap_hit"
            ].fillna(False)
        ]
        .copy()
    )

    multi_tool_summary = (
        build_multi_tool_summary(
            significant_overlaps,
            event_regions,
        )
    )

    ranked_events = (
        build_ranked_splicing_events(
            multi_tool_summary,
            significant_overlaps,
            event_regions,
            overlap_threshold=(
                args.overlap_threshold
            ),
            sashimi_top_n=(
                args.sashimi_top_n
            ),
        )
    )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    paths = {
        "connections": (
            args.outdir
            / (
                f"{args.prefix}"
                ".normalized_connections.tsv"
            )
        ),
        "event_regions": (
            args.outdir
            / (
                f"{args.prefix}"
                ".tool_event_regions.tsv"
            )
        ),
        "all_overlaps": (
            args.outdir
            / (
                f"{args.prefix}"
                ".all_region_overlaps.tsv"
            )
        ),
        "significant_overlaps": (
            args.outdir
            / (
                f"{args.prefix}"
                ".significant_region_overlaps.tsv"
            )
        ),
        "multi_tool_summary": (
            args.outdir
            / (
                f"{args.prefix}"
                ".multi_tool_region_summary.tsv"
            )
        ),
        "ranked_events": (
            args.outdir
            / (
                f"{args.prefix}"
                ".ranked_splicing_events.tsv"
            )
        ),
        "workbook": (
            args.outdir
            / (
                f"{args.prefix}"
                ".splice_region_integration.xlsx"
            )
        ),
    }

    connections.drop(
        columns=[
            "_scope",
            "event_key",
        ],
        errors="ignore",
    ).to_csv(
        paths["connections"],
        sep="\t",
        index=False,
        na_rep="NA",
    )

    event_regions.drop(
        columns=[
            "_scope",
            "event_key",
        ],
        errors="ignore",
    ).to_csv(
        paths["event_regions"],
        sep="\t",
        index=False,
        na_rep="NA",
    )

    all_overlaps.drop(
        columns=[
            "event_key_a",
            "event_key_b",
        ],
        errors="ignore",
    ).to_csv(
        paths["all_overlaps"],
        sep="\t",
        index=False,
        na_rep="NA",
    )

    significant_overlaps.drop(
        columns=[
            "event_key_a",
            "event_key_b",
        ],
        errors="ignore",
    ).to_csv(
        paths[
            "significant_overlaps"
        ],
        sep="\t",
        index=False,
        na_rep="NA",
    )

    multi_tool_summary.to_csv(
        paths[
            "multi_tool_summary"
        ],
        sep="\t",
        index=False,
        na_rep="NA",
    )

    ranked_events.to_csv(
        paths[
            "ranked_events"
        ],
        sep="\t",
        index=False,
        na_rep="NA",
    )

    write_workbook(
        paths["workbook"],
        ranked_events,
        multi_tool_summary,
        significant_overlaps.drop(
            columns=[
                "event_key_a",
                "event_key_b",
            ],
            errors="ignore",
        ),
        all_overlaps.drop(
            columns=[
                "event_key_a",
                "event_key_b",
            ],
            errors="ignore",
        ),
        event_regions,
        connections,
    )

    print(
        f"Connections: "
        f"{len(connections):,}"
    )

    print(
        f"Tool event regions: "
        f"{len(event_regions):,}"
    )

    print(
        f"Overlapping event pairs: "
        f"{len(all_overlaps):,}"
    )

    print(
        f"Significant reciprocal-overlap hits: "
        f"{len(significant_overlaps):,}"
    )

    print(
        f"Multi-tool regions: "
        f"{len(multi_tool_summary):,}"
    )

    print(
        f"Ranked significant regions: "
        f"{len(ranked_events):,}"
    )

    print(
        f"Default sashimi regions: "
        f"{int(ranked_events['default_sashimi_plot'].sum()):,}"
    )

    for label, path in paths.items():
        print(
            f"Wrote {label}: {path}"
        )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(
            main()
        )

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        raise