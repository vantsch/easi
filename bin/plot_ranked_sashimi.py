#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd
import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate ggsashimi plots for top-ranked "
            "EASI splicing events."
        )
    )

    parser.add_argument(
        "--ranked-events",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--bam-manifest",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--gtf",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--workdir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--ggsashimi",
        default="ggsashimi.py",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=30,
    )

    parser.add_argument(
        "--flank",
        type=int,
        default=500,
    )

    parser.add_argument(
        "--min-coverage",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--aggregation",
        default="mean_j",
    )

    parser.add_argument(
        "--output-format",
        default="pdf",
    )

    parser.add_argument(
        "--strand-mode",
        default="NONE",
        choices=[
            "NONE",
            "SENSE",
            "ANTISENSE",
        ],
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.35,
    )

    parser.add_argument(
        "--width",
        type=float,
        default=12.0,
    )

    parser.add_argument(
        "--height",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--annotation-height",
        type=float,
        default=4.5,
    )

    args = parser.parse_args()

    if args.top_n < 1:
        parser.error(
            "--top-n must be at least 1"
        )

    if args.flank < 0:
        parser.error(
            "--flank cannot be negative"
        )

    if args.min_coverage < 0:
        parser.error(
            "--min-coverage cannot be negative"
        )

    if args.width <= 0:
        parser.error(
            "--width must be greater than 0"
        )

    if args.height <= 0:
        parser.error(
            "--height must be greater than 0"
        )

    if args.annotation_height <= 0:
        parser.error(
            "--annotation-height must be greater than 0"
        )

    return args


def require_file(
    path: Path,
    description: str,
) -> Path:
    resolved = (
        path
        .expanduser()
        .resolve()
    )

    if not resolved.is_file():
        raise FileNotFoundError(
            f"{description} not found: {resolved}"
        )

    return resolved


def resolve_executable(
    command: str,
) -> Path:
    supplied = (
        Path(command)
        .expanduser()
    )

    if supplied.is_file():
        return supplied.resolve()

    located = shutil.which(
        command
    )

    if located is not None:
        return Path(
            located
        ).resolve()

    raise FileNotFoundError(
        f"Could not find executable or script: {command}"
    )


def choose_column(
    dataframe: pd.DataFrame,
    candidates: Iterable[str],
    description: str,
) -> str:
    for candidate in candidates:
        if candidate in dataframe.columns:
            return candidate

    raise ValueError(
        f"Could not identify {description} column. "
        f"Expected one of {list(candidates)}; "
        f"found {list(dataframe.columns)}"
    )


def sanitise_filename(
    value: object,
    fallback: str,
) -> str:
    text = str(
        value
    ).strip()

    if (
        not text
        or text.lower() == "nan"
    ):
        text = fallback

    text = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        text,
    )

    text = re.sub(
        r"_+",
        "_",
        text,
    ).strip(
        "._-"
    )

    return text or fallback


def read_ranked_events(
    path: Path,
    top_n: int,
) -> pd.DataFrame:
    events = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    )

    if events.empty:
        raise ValueError(
            f"Ranked-event table is empty: {path}"
        )

    rank_column = choose_column(
        events,
        [
            "rank",
        ],
        "rank",
    )

    chromosome_column = choose_column(
        events,
        [
            "chromosome",
            "chrom",
            "chr",
        ],
        "chromosome",
    )

    start_column = choose_column(
        events,
        [
            "region_start",
            "event_start",
            "start",
        ],
        "event start",
    )

    end_column = choose_column(
        events,
        [
            "region_end",
            "event_end",
            "end",
        ],
        "event end",
    )

    gene_column = next(
        (
            column
            for column in [
                "gene_names",
                "gene_name",
                "gene",
                "genes",
            ]
            if column in events.columns
        ),
        None,
    )

    events = events.copy()

    events["rank"] = pd.to_numeric(
        events[
            rank_column
        ],
        errors="coerce",
    )

    events["chromosome"] = (
        events[
            chromosome_column
        ]
        .astype(str)
        .str.strip()
    )

    events["region_start"] = pd.to_numeric(
        events[
            start_column
        ],
        errors="coerce",
    )

    events["region_end"] = pd.to_numeric(
        events[
            end_column
        ],
        errors="coerce",
    )

    if gene_column is not None:
        events["gene_names"] = (
            events[
                gene_column
            ]
            .fillna(
                "UNANNOTATED"
            )
        )
    else:
        events["gene_names"] = (
            "UNANNOTATED"
        )

    invalid_rows = (
        events[
            [
                "rank",
                "region_start",
                "region_end",
            ]
        ]
        .isna()
        .any(axis=1)
    )

    if invalid_rows.any():
        row_numbers = ", ".join(
            str(index + 2)
            for index in (
                events
                .index[
                    invalid_rows
                ][
                    :10
                ]
            )
        )

        raise ValueError(
            "Non-numeric rank/start/end values "
            f"in ranked-event rows: {row_numbers}"
        )

    events["rank"] = (
        events[
            "rank"
        ]
        .astype(int)
    )

    events["region_start"] = (
        events[
            "region_start"
        ]
        .astype(int)
    )

    events["region_end"] = (
        events[
            "region_end"
        ]
        .astype(int)
    )

    reversed_regions = (
        events[
            "region_start"
        ]
        > events[
            "region_end"
        ]
    )

    if reversed_regions.any():
        original_starts = events.loc[
            reversed_regions,
            "region_start",
        ].copy()

        events.loc[
            reversed_regions,
            "region_start",
        ] = events.loc[
            reversed_regions,
            "region_end",
        ]

        events.loc[
            reversed_regions,
            "region_end",
        ] = original_starts

    selected = (
        events.loc[
            (
                events[
                    "rank"
                ]
                >= 1
            )
            & (
                events[
                    "rank"
                ]
                <= top_n
            )
        ]
        .sort_values(
            [
                "rank",
                "chromosome",
                "region_start",
                "region_end",
            ],
            kind="mergesort",
        )
        .drop_duplicates(
            subset=[
                "rank",
            ],
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    if selected.empty:
        raise ValueError(
            f"No events with ranks 1-{top_n} were found"
        )

    return selected


def read_bam_manifest(
    path: Path,
) -> pd.DataFrame:
    manifest = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    )

    required_columns = [
        "sample_id",
        "bam",
        "condition",
        "bai",
        "label",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in manifest.columns
    ]

    if missing_columns:
        raise ValueError(
            "BAM manifest missing columns: "
            + ", ".join(
                missing_columns
            )
        )

    if manifest.empty:
        raise ValueError(
            f"BAM manifest is empty: {path}"
        )

    for column in required_columns:
        manifest[column] = (
            manifest[
                column
            ]
            .fillna("")
            .astype(str)
            .str.strip()
        )

    for row_number, row in (
        manifest.iterrows()
    ):
        for column in required_columns:
            if not row[column]:
                raise ValueError(
                    f"BAM manifest row {row_number + 2} "
                    f"has an empty {column} value"
                )

    return manifest


def replace_symlink(
    link_path: Path,
    target_path: Path,
) -> None:
    if (
        link_path.exists()
        or link_path.is_symlink()
    ):
        link_path.unlink()

    link_path.symlink_to(
        target_path
    )


def stage_bams(
    manifest: pd.DataFrame,
    workdir: Path,
) -> tuple[
    pd.DataFrame,
    Path,
]:
    bam_directory = (
        workdir
        / "bams"
    )

    bam_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    staged_rows: list[
        dict[str, str]
    ] = []

    ggsashimi_rows: list[
        list[str]
    ] = []

    for index, row in (
        manifest.iterrows()
    ):
        source_bam = require_file(
            Path(
                row[
                    "bam"
                ]
            ),
            (
                "BAM for "
                f"{row['sample_id']}"
            ),
        )

        source_bai = require_file(
            Path(
                row[
                    "bai"
                ]
            ),
            (
                "BAI for "
                f"{row['sample_id']}"
            ),
        )

        sample_token = sanitise_filename(
            row[
                "sample_id"
            ],
            f"sample_{index + 1}",
        )

        staged_bam = (
            bam_directory
            / f"{sample_token}.bam"
        )

        staged_bai = (
            bam_directory
            / f"{sample_token}.bam.bai"
        )

        replace_symlink(
            staged_bam,
            source_bam,
        )

        replace_symlink(
            staged_bai,
            source_bai,
        )

        manifest_bam_path = os.path.relpath(
            staged_bam,
            Path.cwd(),
        )

        manifest_bai_path = os.path.relpath(
            staged_bai,
            Path.cwd(),
        )

        # ggsashimi resolves BAM paths relative to the
        # directory containing ggsashimi_bams.tsv.
        ggsashimi_bam_path = os.path.relpath(
            staged_bam,
            workdir,
        )

        staged_rows.append(
            {
                "sample_id": (
                    row[
                        "sample_id"
                    ]
                ),
                "bam": (
                    manifest_bam_path
                ),
                "condition": (
                    row[
                        "condition"
                    ]
                ),
                "bai": (
                    manifest_bai_path
                ),
                "label": (
                    row[
                        "label"
                    ]
                ),
            }
        )

        ggsashimi_rows.append(
            [
                row[
                    "sample_id"
                ],
                ggsashimi_bam_path,
                row[
                    "condition"
                ],
                row[
                    "label"
                ],
            ]
        )

    staged_manifest = pd.DataFrame(
        staged_rows
    )

    bam_list_path = (
        workdir
        / "ggsashimi_bams.tsv"
    )

    pd.DataFrame(
        ggsashimi_rows
    ).to_csv(
        bam_list_path,
        sep="\t",
        header=False,
        index=False,
    )

    return (
        staged_manifest,
        bam_list_path,
    )


def get_shared_bam_contigs(
    staged_bams: pd.DataFrame,
) -> set[str]:
    contig_sets: list[
        set[str]
    ] = []

    for bam_value in (
        staged_bams[
            "bam"
        ]
    ):
        bam_path = Path(
            str(
                bam_value
            )
        )

        with pysam.AlignmentFile(
            bam_path,
            "rb",
        ) as bam_handle:
            contig_sets.append(
                set(
                    bam_handle.references
                )
            )

    if not contig_sets:
        raise ValueError(
            "No BAM files were available "
            "for contig inspection"
        )

    shared_contigs = (
        contig_sets[
            0
        ]
        .copy()
    )

    for contigs in (
        contig_sets[
            1:
        ]
    ):
        shared_contigs.intersection_update(
            contigs
        )

    if not shared_contigs:
        raise ValueError(
            "The BAM files have no shared contig names"
        )

    return shared_contigs


def contig_candidates(
    chromosome: str,
) -> list[str]:
    chromosome = str(
        chromosome
    ).strip()

    candidates = [
        chromosome,
    ]

    if chromosome.startswith(
        "chr"
    ):
        candidates.append(
            chromosome[
                3:
            ]
        )
    else:
        candidates.append(
            f"chr{chromosome}"
        )

    if chromosome in {
        "chrM",
        "M",
        "MT",
        "chrMT",
    }:
        candidates.extend(
            [
                "chrM",
                "M",
                "MT",
                "chrMT",
            ]
        )

    return list(
        dict.fromkeys(
            candidates
        )
    )


def resolve_contig(
    chromosome: str,
    available_contigs: set[str],
) -> str:
    for candidate in contig_candidates(
        chromosome
    ):
        if candidate in available_contigs:
            return candidate

    raise ValueError(
        f"Could not map chromosome '{chromosome}' "
        "to a BAM contig"
    )


def safe_resolve_contig(
    chromosome: str,
    available_contigs: set[str],
) -> str:
    try:
        return resolve_contig(
            chromosome,
            available_contigs,
        )

    except ValueError:
        return ""


def normalise_gtf_contigs(
    gtf: Path,
    available_contigs: set[str],
    output_path: Path,
) -> Path:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    changed_contigs = 0
    annotation_rows = 0

    with gtf.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as input_handle, output_path.open(
        "w",
        encoding="utf-8",
    ) as output_handle:

        for line in input_handle:
            if (
                line.startswith("#")
                or not line.strip()
            ):
                output_handle.write(
                    line
                )
                continue

            fields = (
                line
                .rstrip(
                    "\n"
                )
                .split(
                    "\t"
                )
            )

            if len(fields) < 9:
                output_handle.write(
                    line
                )
                continue

            original_contig = fields[
                0
            ]

            mapped_contig = safe_resolve_contig(
                original_contig,
                available_contigs,
            )

            if mapped_contig:
                fields[
                    0
                ] = mapped_contig

            if (
                fields[
                    0
                ]
                != original_contig
            ):
                changed_contigs += 1

            output_handle.write(
                "\t".join(
                    fields
                )
                + "\n"
            )

            annotation_rows += 1

    print(
        "Prepared BAM-compatible GTF: "
        f"{annotation_rows:,} rows; "
        f"{changed_contigs:,} contig names changed"
    )

    return output_path


def command_prefix(
    executable: Path,
) -> list[str]:
    if (
        executable
        .suffix
        .lower()
        == ".py"
    ):
        return [
            sys.executable,
            str(
                executable
            ),
        ]

    return [
        str(
            executable
        ),
    ]


def summarise_failure(
    return_code: int,
    stdout: str,
    stderr: str,
    output_exists: bool,
) -> str:
    if (
        return_code == 0
        and not output_exists
    ):
        return (
            "ggsashimi exited successfully, "
            "but the expected plot was absent"
        )

    for text in [
        stderr,
        stdout,
    ]:
        nonempty_lines = [
            line.strip()
            for line in text.splitlines()
            if line.strip()
        ]

        if nonempty_lines:
            return nonempty_lines[
                -1
            ][
                :1000
            ]

    return (
        "ggsashimi exited with "
        f"code {return_code}"
    )


def write_log(
    log_path: Path,
    command: list[str],
    return_code: int,
    stdout: str,
    stderr: str,
) -> None:
    with log_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "COMMAND\n"
        )

        handle.write(
            shlex.join(
                command
            )
            + "\n"
        )

        handle.write(
            f"EXIT_CODE\t{return_code}\n"
        )

        handle.write(
            "STDOUT\n"
        )

        handle.write(
            stdout
        )

        if (
            stdout
            and not stdout.endswith(
                "\n"
            )
        ):
            handle.write(
                "\n"
            )

        handle.write(
            "STDERR\n"
        )

        handle.write(
            stderr
        )

        if (
            stderr
            and not stderr.endswith(
                "\n"
            )
        ):
            handle.write(
                "\n"
            )


def plot_events(
    events: pd.DataFrame,
    *,
    ggsashimi: Path,
    bam_list: Path,
    gtf: Path,
    bam_contigs: set[str],
    outdir: Path,
    flank: int,
    min_coverage: int,
    aggregation: str,
    output_format: str,
    alpha: float,
    width: float,
    height: float,
    annotation_height: float,
    strand_mode: str,
) -> pd.DataFrame:
    results: list[
        dict[str, object]
    ] = []

    for _, event in events.iterrows():
        rank = int(
            event[
                "rank"
            ]
        )

        chromosome = str(
            event[
                "chromosome"
            ]
        ).strip()

        event_start = int(
            event[
                "region_start"
            ]
        )

        event_end = int(
            event[
                "region_end"
            ]
        )

        gene_names = (
            str(
                event[
                    "gene_names"
                ]
            ).strip()
            or "UNANNOTATED"
        )

        plot_start = max(
            1,
            event_start - flank,
        )

        plot_end = (
            event_end
            + flank
        )

        # Limit filename length for events overlapping many genes.
        # Complete gene annotations remain in the ranked-events TSV.
        safe_gene_name = sanitise_filename(
            gene_names,
            "UNANNOTATED",
        )[:50]

        safe_chromosome = sanitise_filename(
            chromosome,
            "contig",
        )[:30]

        output_stem = (
            f"rank_{rank:03d}_"
            f"{safe_gene_name}_"
            f"{safe_chromosome}_"
            f"{event_start}_{event_end}"
        )

        output_prefix = (
            outdir
            / output_stem
        )

        expected_plot = Path(
            f"{output_prefix}.{output_format}"
        )

        log_path = (
            outdir
            / f"{output_stem}.log"
        )

        bam_chromosome = ""
        command: list[str] = []
        stdout = ""
        stderr = ""
        return_code = 1

        try:
            bam_chromosome = resolve_contig(
                chromosome,
                bam_contigs,
            )

            region = (
                f"{bam_chromosome}:"
                f"{plot_start}-"
                f"{plot_end}"
            )

            command = (
                command_prefix(
                    ggsashimi
                )
                + [
                    "-b",
                    str(
                        bam_list
                    ),
                    "-c",
                    region,
                    "-o",
                    str(
                        output_prefix
                    ),
                    "-g",
                    str(
                        gtf
                    ),
                    "--shrink",
                    "-O",
                    "3",
                    "-C",
                    "3",
                    "-L",
                    "4",
                    "-M",
                    str(
                        min_coverage
                    ),
                    "-F",
                    output_format,
                    "--alpha",
                    str(
                        alpha
                    ),
                    "--width",
                    str(
                        width
                    ),
                    "--height",
                    str(
                        height
                    ),
                    "--ann-height",
                    str(
                        annotation_height
                    ),
                    "-s",
                    strand_mode,
                    "-A",
                    aggregation,
                ]
            )

            completed = subprocess.run(
                command,
                check=False,
                capture_output=True,
                text=True,
            )

            return_code = (
                completed.returncode
            )

            stdout = (
                completed.stdout
                or ""
            )

            stderr = (
                completed.stderr
                or ""
            )

        except Exception as error:
            stderr = (
                f"{type(error).__name__}: "
                f"{error}\n"
            )

        output_exists = (
            expected_plot.is_file()
        )

        success = (
            return_code == 0
            and output_exists
        )

        status = (
            "SUCCESS"
            if success
            else "FAILED"
        )

        error_message = ""

        if not success:
            error_message = summarise_failure(
                return_code,
                stdout,
                stderr,
                output_exists,
            )

        write_log(
            log_path,
            command,
            return_code,
            stdout,
            stderr,
        )

        results.append(
            {
                "rank": (
                    rank
                ),
                "gene_names": (
                    gene_names
                ),
                "chromosome": (
                    chromosome
                ),
                "bam_chromosome": (
                    bam_chromosome
                ),
                "event_start": (
                    event_start
                ),
                "event_end": (
                    event_end
                ),
                "plot_start": (
                    plot_start
                ),
                "plot_end": (
                    plot_end
                ),
                "status": (
                    status
                ),
                "exit_code": (
                    return_code
                ),
                "plot_file": (
                    str(
                        expected_plot
                    )
                    if output_exists
                    else ""
                ),
                "log_file": (
                    str(
                        log_path
                    )
                ),
                "error": (
                    error_message
                ),
            }
        )

        message = (
            f"[{status}] rank {rank}: "
            f"{gene_names} "
            f"{chromosome}:"
            f"{plot_start}-"
            f"{plot_end}"
        )

        print(
            message,
            file=(
                sys.stdout
                if success
                else sys.stderr
            ),
        )

    return pd.DataFrame(
        results
    )


def main() -> int:
    args = parse_args()

    ranked_events_path = require_file(
        args.ranked_events,
        "Ranked-event table",
    )

    bam_manifest_path = require_file(
        args.bam_manifest,
        "BAM manifest",
    )

    gtf_path = require_file(
        args.gtf,
        "GTF annotation",
    )

    ggsashimi_path = resolve_executable(
        args.ggsashimi
    )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    selected_events = read_ranked_events(
        ranked_events_path,
        args.top_n,
    )

    bam_manifest = read_bam_manifest(
        bam_manifest_path
    )

    (
        staged_bams,
        ggsashimi_bam_list,
    ) = stage_bams(
        bam_manifest,
        args.workdir,
    )

    bam_contigs = get_shared_bam_contigs(
        staged_bams
    )

    selected_events[
        "bam_chromosome"
    ] = selected_events[
        "chromosome"
    ].map(
        lambda chromosome: safe_resolve_contig(
            chromosome,
            bam_contigs,
        )
    )

    compatible_gtf = normalise_gtf_contigs(
        gtf_path,
        bam_contigs,
        (
            args.workdir
            / "annotation.bam_contigs.gtf"
        ),
    )

    selected_events_path = (
        args.outdir
        / "selected_ranked_events.tsv"
    )

    staged_manifest_path = (
        args.outdir
        / "sashimi_bam_manifest.tsv"
    )

    sashimi_manifest_path = (
        args.outdir
        / "sashimi_manifest.tsv"
    )

    selected_events.to_csv(
        selected_events_path,
        sep="\t",
        index=False,
    )

    staged_bams.to_csv(
        staged_manifest_path,
        sep="\t",
        index=False,
    )

    sashimi_manifest = plot_events(
        selected_events,
        ggsashimi=(
            ggsashimi_path
        ),
        bam_list=(
            ggsashimi_bam_list
        ),
        gtf=(
            compatible_gtf
        ),
        bam_contigs=(
            bam_contigs
        ),
        outdir=(
            args.outdir
        ),
        flank=(
            args.flank
        ),
        min_coverage=(
            args.min_coverage
        ),
        aggregation=(
            args.aggregation
        ),
        output_format=(
            args.output_format
        ),
        alpha=(
            args.alpha
        ),
        width=(
            args.width
        ),
        height=(
            args.height
        ),
        annotation_height=(
            args.annotation_height
        ),
        strand_mode=(
            args.strand_mode
        ),
    )

    sashimi_manifest.to_csv(
        sashimi_manifest_path,
        sep="\t",
        index=False,
    )

    successful = int(
        (
            sashimi_manifest[
                "status"
            ]
            == "SUCCESS"
        ).sum()
    )

    failed = int(
        (
            sashimi_manifest[
                "status"
            ]
            == "FAILED"
        ).sum()
    )

    print(
        "Selected ranked events: "
        f"{len(selected_events)}"
    )

    print(
        "Successful sashimi plots: "
        f"{successful}"
    )

    print(
        "Failed sashimi plots: "
        f"{failed}"
    )

    print(
        "Wrote manifest: "
        f"{sashimi_manifest_path}"
    )

    return (
        0
        if successful > 0
        else 1
    )


if __name__ == "__main__":
    try:
        raise SystemExit(
            main()
        )

    except Exception as error:
        print(
            "ERROR: "
            f"{type(error).__name__}: "
            f"{error}",
            file=sys.stderr,
        )

        raise SystemExit(
            1
        )