#!/usr/bin/env python3

"""
Integrate an ordinary DESeq2 results table with standardized MAJIQ, rMATS,
LeafCutter and IsoformSwitchAnalyzeR gene-summary tables.

DESeq2 input:
    The ordinary unfiltered DESeq2 results table, for example:

        UHR_vs_HBR.deseq2.results.tsv

    The DESeq2 table is converted internally to the same gene-level structure
    used by the standardized splicing summaries.

Splicing inputs:
    MAJIQ, rMATS, LeafCutter and ISA must each contain one standardized
    gene-level row per gene with the following columns:

        gene_id
        gene_name
        method
        detected
        significant
        n_events
        n_significant_events
        max_abs_effect
        effect_direction
        event_types

DESeq2 is used as the integration backbone because its unfiltered results
table should contain the most complete tested gene universe. MAJIQ, rMATS,
LeafCutter and ISA are left-joined onto the DESeq2 genes.

Genes found only by a splicing method and absent from the DESeq2 results are
therefore not added as new rows.

Outputs:
    <prefix>.gene_integration_short.tsv
    <prefix>.gene_integration_long.tsv
    <prefix>.gene_ranking_parameters.tsv
    <prefix>.gene_ranking_parameters.json

Ranking:
    expression_score =
        category_weight_expression
        * method_weight_deseq2
        * deseq2_significant

    splicing_score =
        category_weight_splicing
        * weighted mean of significant MAJIQ, rMATS, LeafCutter and ISA calls

    agreement_bonus =
        min(
            max(n_methods_significant - 1, 0)
            * per_additional_method,
            maximum
        )

    gene_rank_score =
        expression_score
        + splicing_score
        + agreement_bonus

Rows are sorted only after the ranking metrics have been calculated.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


DEFAULT_CONFIG: dict[str, Any] = {
    "thresholds": {
        "deseq2": {
            "padj": 0.05,
            "abs_log2fc": 0.0,
        },
    },
    "ranking": {
        "method_weights": {
            "deseq2": 1.0,
            "majiq": 1.0,
            "rmats": 1.0,
            "leafcutter": 1.0,
            "isa": 1.0,
        },
        "category_weights": {
            "expression": 1.0,
            "splicing": 1.0,
        },
        "agreement_bonus": {
            "enabled": True,
            "per_additional_method": 0.20,
            "maximum": 0.60,
        },
    },
}


METHODS = [
    "deseq2",
    "majiq",
    "rmats",
    "leafcutter",
    "isa",
]

SPLICING_METHODS = [
    "majiq",
    "rmats",
    "leafcutter",
    "isa",
]

DISPLAY_NAMES = {
    "deseq2": "DESeq2",
    "majiq": "MAJIQ",
    "rmats": "rMATS",
    "leafcutter": "LeafCutter",
    "isa": "ISA",
}

EXPECTED_METHOD_NAMES = {
    "majiq": {
        "majiq",
    },
    "rmats": {
        "rmats",
    },
    "leafcutter": {
        "leafcutter",
    },
    "isa": {
        "isa",
        "isoformswitchanalyzer",
        "isoform_switch_analyzer",
    },
}

REQUIRED_COLUMNS = [
    "gene_id",
    "gene_name",
    "method",
    "detected",
    "significant",
    "n_events",
    "n_significant_events",
    "max_abs_effect",
    "effect_direction",
    "event_types",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Integrate an ordinary DESeq2 results table with standardized "
            "EASI splicing gene-summary tables, using DESeq2 as the "
            "left-join backbone."
        )
    )

    parser.add_argument(
        "--deseq2",
        required=True,
        type=Path,
        help=(
            "Ordinary unfiltered DESeq2 results table, for example "
            "'UHR_vs_HBR.deseq2.results.tsv'."
        ),
    )

    parser.add_argument(
        "--majiq",
        required=True,
        type=Path,
        help="Standardized MAJIQ gene-summary table.",
    )

    parser.add_argument(
        "--rmats",
        required=True,
        type=Path,
        help="Standardized rMATS gene-summary table.",
    )

    parser.add_argument(
        "--leafcutter",
        required=True,
        type=Path,
        help="Standardized LeafCutter gene-summary table.",
    )

    parser.add_argument(
        "--isa",
        required=True,
        type=Path,
        help="Standardized ISA gene-summary table.",
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help=(
            "Optional YAML file containing DESeq2 thresholds and "
            "gene-ranking settings."
        ),
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Output prefix, for example UHR_vs_HBR.",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Output directory.",
    )

    return parser.parse_args()


def normalise_column_name(value: str) -> str:
    value = str(value).strip().lower()

    value = re.sub(
        r"[^a-z0-9]+",
        "_",
        value,
    )

    return value.strip("_")


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Input file does not exist: {path}"
        )

    if path.suffix.lower() in {
        ".xlsx",
        ".xls",
    }:
        dataframe = pd.read_excel(
            path,
            dtype=str,
        )
    else:
        dataframe = pd.read_csv(
            path,
            sep=None,
            engine="python",
            dtype=str,
        )

    dataframe.columns = [
        normalise_column_name(column)
        for column in dataframe.columns
    ]

    return dataframe


def merge_config(
    default: dict[str, Any],
    custom: dict[str, Any],
) -> dict[str, Any]:
    result = deepcopy(
        default
    )

    for key, value in custom.items():
        if (
            key in result
            and isinstance(result[key], dict)
            and isinstance(value, dict)
        ):
            result[key] = merge_config(
                result[key],
                value,
            )
        else:
            result[key] = value

    return result


def load_config(
    path: Path | None,
) -> dict[str, Any]:
    if path is None:
        return deepcopy(
            DEFAULT_CONFIG
        )

    if not path.exists():
        raise FileNotFoundError(
            f"Configuration file does not exist: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        custom = yaml.safe_load(
            handle
        ) or {}

    if not isinstance(
        custom,
        dict,
    ):
        raise ValueError(
            f"Configuration root must be a YAML mapping: {path}"
        )

    return merge_config(
        DEFAULT_CONFIG,
        custom,
    )


def validate_config(
    config: dict[str, Any],
) -> None:
    try:
        padj_threshold = float(
            config[
                "thresholds"
            ][
                "deseq2"
            ][
                "padj"
            ]
        )

        abs_log2fc_threshold = float(
            config[
                "thresholds"
            ][
                "deseq2"
            ][
                "abs_log2fc"
            ]
        )

    except (
        KeyError,
        TypeError,
        ValueError,
    ) as error:
        raise ValueError(
            "Invalid DESeq2 threshold configuration. Expected:\n"
            "thresholds:\n"
            "  deseq2:\n"
            "    padj: 0.05\n"
            "    abs_log2fc: 0.0"
        ) from error

    if not 0 <= padj_threshold <= 1:
        raise ValueError(
            "thresholds.deseq2.padj must be between 0 and 1."
        )

    if abs_log2fc_threshold < 0:
        raise ValueError(
            "thresholds.deseq2.abs_log2fc must be greater than "
            "or equal to 0."
        )

    try:
        method_weights = config[
            "ranking"
        ][
            "method_weights"
        ]

        category_weights = config[
            "ranking"
        ][
            "category_weights"
        ]

        bonus_config = config[
            "ranking"
        ][
            "agreement_bonus"
        ]

    except (
        KeyError,
        TypeError,
    ) as error:
        raise ValueError(
            "The ranking configuration is incomplete."
        ) from error

    for method in METHODS:
        if method not in method_weights:
            raise ValueError(
                "Missing ranking method weight for "
                f"'{method}'."
            )

        try:
            weight = float(
                method_weights[
                    method
                ]
            )

        except (
            TypeError,
            ValueError,
        ) as error:
            raise ValueError(
                f"Ranking weight for '{method}' must be numeric."
            ) from error

        if weight < 0:
            raise ValueError(
                f"Ranking weight for '{method}' cannot be negative."
            )

    for category in [
        "expression",
        "splicing",
    ]:
        if category not in category_weights:
            raise ValueError(
                "Missing ranking category weight for "
                f"'{category}'."
            )

        try:
            weight = float(
                category_weights[
                    category
                ]
            )

        except (
            TypeError,
            ValueError,
        ) as error:
            raise ValueError(
                f"Category weight for '{category}' must be numeric."
            ) from error

        if weight < 0:
            raise ValueError(
                f"Category weight for '{category}' cannot be negative."
            )

    for parameter in [
        "enabled",
        "per_additional_method",
        "maximum",
    ]:
        if parameter not in bonus_config:
            raise ValueError(
                "Missing agreement-bonus parameter "
                f"'{parameter}'."
            )

    try:
        per_additional_method = float(
            bonus_config[
                "per_additional_method"
            ]
        )

        maximum = float(
            bonus_config[
                "maximum"
            ]
        )

    except (
        TypeError,
        ValueError,
    ) as error:
        raise ValueError(
            "Agreement-bonus values must be numeric."
        ) from error

    if per_additional_method < 0:
        raise ValueError(
            "ranking.agreement_bonus.per_additional_method "
            "cannot be negative."
        )

    if maximum < 0:
        raise ValueError(
            "ranking.agreement_bonus.maximum cannot be negative."
        )


def clean_gene_id(
    series: pd.Series,
) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.replace(
            r"\.\d+$",
            "",
            regex=True,
        )
        .replace(
            {
                "": pd.NA,
                "NA": pd.NA,
                "nan": pd.NA,
                "None": pd.NA,
            }
        )
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
                "nan": pd.NA,
                "None": pd.NA,
            }
        )
    )


def parse_boolean_value(
    value: Any,
) -> Any:
    if pd.isna(
        value
    ):
        return pd.NA

    text = normalise_column_name(
        str(value)
    )

    true_values = {
        "true",
        "t",
        "yes",
        "y",
        "1",
        "significant",
        "detected",
    }

    false_values = {
        "false",
        "f",
        "no",
        "n",
        "0",
        "not_significant",
        "nonsignificant",
        "not_detected",
        "undetected",
    }

    if text in true_values:
        return True

    if text in false_values:
        return False

    return pd.NA


def to_boolean(
    series: pd.Series,
) -> pd.Series:
    return (
        series.map(
            parse_boolean_value
        )
        .astype("boolean")
    )


def to_integer(
    series: pd.Series,
) -> pd.Series:
    return (
        pd.to_numeric(
            series,
            errors="coerce",
        )
        .round()
        .astype("Int64")
    )


def to_float(
    series: pd.Series,
) -> pd.Series:
    return pd.to_numeric(
        series,
        errors="coerce",
    )


def first_nonempty(
    series: pd.Series,
) -> Any:
    for value in series:
        if pd.notna(
            value
        ):
            text = str(
                value
            ).strip()

            if text:
                return text

    return pd.NA


def combine_unique_text(
    series: pd.Series,
    separator: str = ";",
) -> Any:
    values: list[str] = []

    for value in series:
        if pd.isna(
            value
        ):
            continue

        for item in str(
            value
        ).split(
            separator
        ):
            item = item.strip()

            if (
                item
                and item not in values
            ):
                values.append(
                    item
                )

    if not values:
        return pd.NA

    return separator.join(
        values
    )


def select_effect_direction(
    dataframe: pd.DataFrame,
) -> Any:
    significant_values = dataframe.loc[
        dataframe[
            "significant"
        ].fillna(
            False
        ),
        "effect_direction",
    ]

    value = first_nonempty(
        significant_values
    )

    if pd.notna(
        value
    ):
        return value

    return first_nonempty(
        dataframe[
            "effect_direction"
        ]
    )


def find_column(
    dataframe: pd.DataFrame,
    aliases: list[str],
    required: bool = True,
) -> str | None:
    """
    Find the first matching normalized column name.
    """

    normalized_aliases = [
        normalise_column_name(
            alias
        )
        for alias in aliases
    ]

    for alias in normalized_aliases:
        if alias in dataframe.columns:
            return alias

    if required:
        raise ValueError(
            "Could not find any of the expected columns: "
            f"{', '.join(normalized_aliases)}. "
            "Available columns: "
            f"{', '.join(map(str, dataframe.columns))}"
        )

    return None


def collapse_duplicate_genes(
    dataframe: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    """
    Defensively collapse duplicate gene IDs to prevent many-to-many joins.

    For duplicate rows:
      - detected and significant are combined with any()
      - event counts are summed
      - max_abs_effect is maximized
      - event_types are combined
      - effect_direction is taken from a significant row where possible
    """

    duplicate_mask = dataframe[
        "gene_id"
    ].duplicated(
        keep=False
    )

    if not duplicate_mask.any():
        return dataframe

    duplicate_count = int(
        dataframe.loc[
            duplicate_mask,
            "gene_id",
        ].nunique()
    )

    print(
        f"WARNING: {DISPLAY_NAMES[method]} contains "
        f"{duplicate_count} duplicated gene IDs. "
        "Collapsing them before integration.",
        file=sys.stderr,
    )

    rows: list[dict[str, Any]] = []

    for gene_id, gene_data in dataframe.groupby(
        "gene_id",
        sort=False,
        dropna=False,
    ):
        max_abs_effect = gene_data[
            "max_abs_effect"
        ].max(
            skipna=True
        )

        if pd.isna(
            max_abs_effect
        ):
            max_abs_effect = np.nan

        rows.append(
            {
                "gene_id": gene_id,
                "gene_name": first_nonempty(
                    gene_data[
                        "gene_name"
                    ]
                ),
                "method": first_nonempty(
                    gene_data[
                        "method"
                    ]
                ),
                "detected": bool(
                    gene_data[
                        "detected"
                    ]
                    .fillna(False)
                    .any()
                ),
                "significant": bool(
                    gene_data[
                        "significant"
                    ]
                    .fillna(False)
                    .any()
                ),
                "n_events": int(
                    gene_data[
                        "n_events"
                    ]
                    .fillna(0)
                    .sum()
                ),
                "n_significant_events": int(
                    gene_data[
                        "n_significant_events"
                    ]
                    .fillna(0)
                    .sum()
                ),
                "max_abs_effect": max_abs_effect,
                "effect_direction": select_effect_direction(
                    gene_data
                ),
                "event_types": combine_unique_text(
                    gene_data[
                        "event_types"
                    ]
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


def load_deseq2_results(
    path: Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Load an ordinary unfiltered DESeq2 results table and convert it to the
    common internal gene-level schema.

    Required information:
      - gene ID
      - base mean
      - log2 fold change
      - adjusted p-value

    Gene name is optional. If absent, it can later be filled from a splicing
    summary during integration.
    """

    dataframe = read_table(
        path
    )

    gene_id_column = find_column(
        dataframe,
        [
            "gene_id",
            "gene",
            "ensembl_gene_id",
            "ensembl",
            "geneid",
        ],
    )

    gene_name_column = find_column(
        dataframe,
        [
            "gene_name",
            "gene_symbol",
            "symbol",
            "external_gene_name",
            "genename",
        ],
        required=False,
    )

    base_mean_column = find_column(
        dataframe,
        [
            "baseMean",
            "base_mean",
            "basemean",
        ],
    )

    log2fc_column = find_column(
        dataframe,
        [
            "log2FoldChange",
            "log2_fold_change",
            "log2fc",
            "logfc",
        ],
    )

    padj_column = find_column(
        dataframe,
        [
            "padj",
            "adjusted_p_value",
            "adjusted_pvalue",
            "fdr",
            "qvalue",
            "q_value",
        ],
    )

    gene_id = clean_gene_id(
        dataframe[
            gene_id_column
        ]
    )

    if gene_name_column is not None:
        gene_name = clean_text(
            dataframe[
                gene_name_column
            ]
        )
    else:
        gene_name = pd.Series(
            pd.NA,
            index=dataframe.index,
            dtype="string",
        )

    base_mean = to_float(
        dataframe[
            base_mean_column
        ]
    )

    log2fc = to_float(
        dataframe[
            log2fc_column
        ]
    )

    padj = to_float(
        dataframe[
            padj_column
        ]
    )

    thresholds = config[
        "thresholds"
    ][
        "deseq2"
    ]

    padj_threshold = float(
        thresholds[
            "padj"
        ]
    )

    abs_log2fc_threshold = float(
        thresholds[
            "abs_log2fc"
        ]
    )

    detected = (
        base_mean.notna()
        & base_mean.gt(0)
    )

    significant = (
        padj.notna()
        & padj.lt(
            padj_threshold
        )
        & log2fc.notna()
        & log2fc.abs().ge(
            abs_log2fc_threshold
        )
    )

    result = pd.DataFrame(
        {
            "gene_id": gene_id,
            "gene_name": gene_name,
            "method": "deseq2",
            "detected": detected,
            "significant": significant,
            "log2fc": log2fc,
        }
    )

    result = result[
        result[
            "gene_id"
        ].notna()
    ].copy()

    result[
        "detected"
    ] = (
        result[
            "detected"
        ]
        .fillna(False)
        .astype(bool)
    )

    result[
        "significant"
    ] = (
        result[
            "significant"
        ]
        .fillna(False)
        .astype(bool)
    )

    if result[
        "gene_id"
    ].duplicated().any():
        duplicated = result.loc[
            result[
                "gene_id"
            ].duplicated(
                keep=False
            ),
            "gene_id",
        ].unique()

        raise ValueError(
            "DESeq2 results contain duplicated gene IDs after ID normalization: "
            + ", ".join(
                map(
                    str,
                    duplicated[:10],
                )
            )
        )

    return result.reset_index(
        drop=True
    )


def validate_required_columns(
    dataframe: pd.DataFrame,
    method: str,
    path: Path,
) -> None:
    missing = [
        column
        for column in REQUIRED_COLUMNS
        if column not in dataframe.columns
    ]

    if missing:
        available = ", ".join(
            map(
                str,
                dataframe.columns,
            )
        )

        raise ValueError(
            f"{DISPLAY_NAMES[method]} gene summary is missing "
            f"required columns: {', '.join(missing)}. "
            f"File: {path}. "
            f"Available columns: {available}"
        )


def validate_method_column(
    dataframe: pd.DataFrame,
    method: str,
    path: Path,
) -> None:
    observed = {
        normalise_column_name(
            value
        )
        for value in dataframe[
            "method"
        ].dropna().unique()
        if str(
            value
        ).strip()
    }

    if not observed:
        return

    expected = EXPECTED_METHOD_NAMES[
        method
    ]

    if not observed.issubset(
        expected
    ):
        raise ValueError(
            f"Unexpected values in the method column for "
            f"{DISPLAY_NAMES[method]}: {sorted(observed)}. "
            f"Expected one of: {sorted(expected)}. "
            f"File: {path}"
        )


def load_gene_summary(
    path: Path,
    method: str,
) -> pd.DataFrame:
    dataframe = read_table(
        path
    )

    validate_required_columns(
        dataframe,
        method,
        path,
    )

    dataframe = dataframe[
        REQUIRED_COLUMNS
    ].copy()

    dataframe[
        "gene_id"
    ] = clean_gene_id(
        dataframe[
            "gene_id"
        ]
    )

    dataframe[
        "gene_name"
    ] = clean_text(
        dataframe[
            "gene_name"
        ]
    )

    dataframe[
        "method"
    ] = clean_text(
        dataframe[
            "method"
        ]
    )

    dataframe[
        "detected"
    ] = to_boolean(
        dataframe[
            "detected"
        ]
    )

    dataframe[
        "significant"
    ] = to_boolean(
        dataframe[
            "significant"
        ]
    )

    dataframe[
        "n_events"
    ] = to_integer(
        dataframe[
            "n_events"
        ]
    )

    dataframe[
        "n_significant_events"
    ] = to_integer(
        dataframe[
            "n_significant_events"
        ]
    )

    dataframe[
        "max_abs_effect"
    ] = to_float(
        dataframe[
            "max_abs_effect"
        ]
    )

    dataframe[
        "effect_direction"
    ] = clean_text(
        dataframe[
            "effect_direction"
        ]
    )

    dataframe[
        "event_types"
    ] = clean_text(
        dataframe[
            "event_types"
        ]
    )

    dataframe = dataframe[
        dataframe[
            "gene_id"
        ].notna()
    ].copy()

    validate_method_column(
        dataframe,
        method,
        path,
    )

    dataframe[
        "detected"
    ] = (
        dataframe[
            "detected"
        ]
        .fillna(False)
        .astype(bool)
    )

    dataframe[
        "significant"
    ] = (
        dataframe[
            "significant"
        ]
        .fillna(False)
        .astype(bool)
    )

    dataframe[
        "n_events"
    ] = (
        dataframe[
            "n_events"
        ]
        .fillna(0)
        .astype(int)
    )

    dataframe[
        "n_significant_events"
    ] = (
        dataframe[
            "n_significant_events"
        ]
        .fillna(0)
        .astype(int)
    )

    dataframe = collapse_duplicate_genes(
        dataframe,
        method,
    )

    return dataframe.reset_index(
        drop=True
    )


def rename_method_columns(
    dataframe: pd.DataFrame,
    method: str,
    keep_gene_name: bool = False,
) -> pd.DataFrame:
    rename_map = {
        "method": f"{method}_method",
        "detected": f"{method}_detected",
        "significant": f"{method}_significant",
        "n_events": f"{method}_n_events",
        "n_significant_events": f"{method}_n_significant_events",
        "max_abs_effect": f"{method}_max_abs_effect",
        "effect_direction": f"{method}_effect_direction",
        "log2fc": f"{method}_log2fc",
        "event_types": f"{method}_event_types",
    }

    if not keep_gene_name:
        rename_map[
            "gene_name"
        ] = f"{method}_gene_name"

    return dataframe.rename(
        columns=rename_map
    )


def integrate_methods(
    method_tables: dict[str, pd.DataFrame],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Use DESeq2 as the fixed gene universe and left-join all standardized
    splicing summaries onto it.
    """

    integrated = rename_method_columns(
        method_tables[
            "deseq2"
        ].copy(),
        "deseq2",
        keep_gene_name=True,
    )

    for method in SPLICING_METHODS:
        method_table = rename_method_columns(
            method_tables[
                method
            ].copy(),
            method,
            keep_gene_name=False,
        )

        integrated = integrated.merge(
            method_table,
            on="gene_id",
            how="left",
            validate="one_to_one",
            sort=False,
        )

        method_gene_name = (
            f"{method}_gene_name"
        )

        integrated[
            "gene_name"
        ] = (
            integrated[
                "gene_name"
            ]
            .fillna(
                integrated[
                    method_gene_name
                ]
            )
        )

        integrated = integrated.drop(
            columns=[
                method_gene_name
            ]
        )

    for method in METHODS:
        detected_column = (
            f"{method}_detected"
        )

        significant_column = (
            f"{method}_significant"
        )

        integrated[
            detected_column
        ] = (
            integrated[
                detected_column
            ]
            .fillna(False)
            .astype(bool)
        )

        integrated[
            significant_column
        ] = (
            integrated[
                significant_column
            ]
            .fillna(False)
            .astype(bool)
        )

    for method in SPLICING_METHODS:
        n_events_column = (
            f"{method}_n_events"
        )

        n_significant_column = (
            f"{method}_n_significant_events"
        )

        integrated[
            n_events_column
        ] = (
            integrated[
                n_events_column
            ]
            .fillna(0)
            .astype(int)
        )

        integrated[
            n_significant_column
        ] = (
            integrated[
                n_significant_column
            ]
            .fillna(0)
            .astype(int)
        )

    integrated[
        "expression_significant"
    ] = integrated[
        "deseq2_significant"
    ]

    integrated[
        "splicing_significant"
    ] = integrated[
        [
            f"{method}_significant"
            for method in SPLICING_METHODS
        ]
    ].any(
        axis=1
    )

    integrated[
        "n_methods_detected"
    ] = integrated[
        [
            f"{method}_detected"
            for method in METHODS
        ]
    ].sum(
        axis=1
    ).astype(
        int
    )

    integrated[
        "n_methods_significant"
    ] = integrated[
        [
            f"{method}_significant"
            for method in METHODS
        ]
    ].sum(
        axis=1
    ).astype(
        int
    )

    integrated[
        "significant_methods"
    ] = integrated.apply(
        lambda row: ";".join(
            DISPLAY_NAMES[
                method
            ]
            for method in METHODS
            if row[
                f"{method}_significant"
            ]
        ),
        axis=1,
    )

    method_weights = config[
        "ranking"
    ][
        "method_weights"
    ]

    category_weights = config[
        "ranking"
    ][
        "category_weights"
    ]

    for method in METHODS:
        integrated[
            f"{method}_rank_component"
        ] = (
            integrated[
                f"{method}_significant"
            ]
            .astype(float)
            .mul(
                float(
                    method_weights[
                        method
                    ]
                )
            )
        )

    integrated[
        "expression_category_score"
    ] = (
        integrated[
            "deseq2_rank_component"
        ]
        * float(
            category_weights[
                "expression"
            ]
        )
    )

    splicing_weight_sum = sum(
        float(
            method_weights[
                method
            ]
        )
        for method in SPLICING_METHODS
    )

    if splicing_weight_sum > 0:
        splicing_components = integrated[
            [
                f"{method}_rank_component"
                for method in SPLICING_METHODS
            ]
        ].sum(
            axis=1
        )

        integrated[
            "splicing_category_score"
        ] = (
            splicing_components
            / splicing_weight_sum
            * float(
                category_weights[
                    "splicing"
                ]
            )
        )

    else:
        integrated[
            "splicing_category_score"
        ] = 0.0

    bonus_config = config[
        "ranking"
    ][
        "agreement_bonus"
    ]

    if bool(
        bonus_config[
            "enabled"
        ]
    ):
        integrated[
            "agreement_bonus"
        ] = (
            integrated[
                "n_methods_significant"
            ]
            .sub(1)
            .clip(
                lower=0
            )
            .mul(
                float(
                    bonus_config[
                        "per_additional_method"
                    ]
                )
            )
            .clip(
                upper=float(
                    bonus_config[
                        "maximum"
                    ]
                )
            )
        )

    else:
        integrated[
            "agreement_bonus"
        ] = 0.0

    integrated[
        "gene_rank_score"
    ] = (
        integrated[
            "expression_category_score"
        ]
        + integrated[
            "splicing_category_score"
        ]
        + integrated[
            "agreement_bonus"
        ]
    )

    integrated = integrated.sort_values(
        by=[
            "gene_rank_score",
            "n_methods_significant",
            "splicing_category_score",
            "expression_category_score",
            "gene_name",
            "gene_id",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            True,
            True,
        ],
        na_position="last",
        kind="mergesort",
    ).reset_index(
        drop=True
    )

    integrated.insert(
        0,
        "gene_rank",
        np.arange(
            1,
            len(
                integrated
            ) + 1,
        ),
    )

    return integrated


def make_short_table(
    integrated: pd.DataFrame,
) -> pd.DataFrame:
    short = integrated[
        [
            "gene_rank",
            "gene_id",
            "gene_name",
            "gene_rank_score",
            "deseq2_significant",
            "majiq_significant",
            "rmats_significant",
            "leafcutter_significant",
            "isa_significant",
            "n_methods_significant",
            "significant_methods",
        ]
    ].copy()

    short = short.rename(
        columns={
            "gene_name": "gene_symbol",
            "deseq2_significant": "deseq2",
            "majiq_significant": "majiq",
            "rmats_significant": "rmats",
            "leafcutter_significant": "leafcutter",
            "isa_significant": "isa",
        }
    )

    for method in METHODS:
        short[
            method
        ] = short[
            method
        ].map(
            {
                True: "yes",
                False: "no",
            }
        )

    return short


def make_long_table(
    integrated: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create the comprehensive gene-level table used by easiR.

    Internal helper columns used to construct the ranking, such as
    *_method and *_rank_component, are deliberately not exported.
    """

    columns = [
        "gene_rank",
        "gene_id",
        "gene_name",
        "gene_rank_score",

        "deseq2_detected",
        "deseq2_significant",
        "deseq2_log2fc",

        "majiq_detected",
        "majiq_significant",
        "majiq_n_events",
        "majiq_n_significant_events",
        "majiq_max_abs_effect",
        "majiq_effect_direction",
        "majiq_event_types",

        "rmats_detected",
        "rmats_significant",
        "rmats_n_events",
        "rmats_n_significant_events",
        "rmats_max_abs_effect",
        "rmats_effect_direction",
        "rmats_event_types",

        "leafcutter_detected",
        "leafcutter_significant",
        "leafcutter_n_events",
        "leafcutter_n_significant_events",
        "leafcutter_max_abs_effect",
        "leafcutter_effect_direction",
        "leafcutter_event_types",

        "isa_detected",
        "isa_significant",
        "isa_n_events",
        "isa_n_significant_events",
        "isa_max_abs_effect",
        "isa_effect_direction",
        "isa_event_types",

        "expression_significant",
        "splicing_significant",
        "n_methods_detected",
        "n_methods_significant",
        "significant_methods",

        "expression_category_score",
        "splicing_category_score",
        "agreement_bonus",
    ]

    return (
        integrated[
            columns
        ]
        .rename(
            columns={
                "gene_name": "gene_symbol",
            }
        )
        .copy()
    )


def flatten_config(
    value: Any,
    prefix: str = "",
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    if isinstance(
        value,
        dict,
    ):
        for key, child in value.items():
            child_prefix = (
                f"{prefix}.{key}"
                if prefix
                else key
            )

            rows.extend(
                flatten_config(
                    child,
                    child_prefix,
                )
            )

    else:
        rows.append(
            {
                "parameter": prefix,
                "value": value,
            }
        )

    return rows


def write_gene_integration_tsv(
    dataframe: pd.DataFrame,
    path: Path,
    prefix: str,
) -> None:
    metadata = [
        "# EASI gene integration ranking",
        (
            "# expression_category_score = "
            "category_weight_expression * method_weight_deseq2 "
            "* deseq2_significant"
        ),
        (
            "# splicing_category_score = "
            "category_weight_splicing * weighted mean of significant calls "
            "from MAJIQ, rMATS, LeafCutter and ISA"
        ),
        (
            "# agreement_bonus = "
            "min(max(n_methods_significant - 1, 0) "
            "* per_additional_method, maximum)"
        ),
        (
            "# gene_rank_score = "
            "expression_category_score + splicing_category_score "
            "+ agreement_bonus"
        ),
        (
            "# rank_order = gene_rank_score desc; "
            "n_methods_significant desc; "
            "splicing_category_score desc; "
            "expression_category_score desc; "
            "gene_name asc; gene_id asc"
        ),
        (
            "# Ranking parameter values are reported in "
            f"{prefix}.gene_ranking_parameters.tsv"
        ),
    ]

    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for line in metadata:
            handle.write(
                line + "\n"
            )

        dataframe.to_csv(
            handle,
            sep="\t",
            index=False,
            na_rep="NA",
        )


def main() -> int:
    args = parse_arguments()

    config = load_config(
        args.config
    )

    validate_config(
        config
    )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    method_tables = {
        "deseq2": load_deseq2_results(
            args.deseq2,
            config,
        ),

        "majiq": load_gene_summary(
            args.majiq,
            "majiq",
        ),

        "rmats": load_gene_summary(
            args.rmats,
            "rmats",
        ),

        "leafcutter": load_gene_summary(
            args.leafcutter,
            "leafcutter",
        ),

        "isa": load_gene_summary(
            args.isa,
            "isa",
        ),
    }

    integrated = integrate_methods(
        method_tables,
        config,
    )

    short_table = make_short_table(
        integrated
    )

    long_table = make_long_table(
        integrated
    )

    short_path = (
        args.outdir
        / (
            f"{args.prefix}"
            ".gene_integration_short.tsv"
        )
    )

    long_path = (
        args.outdir
        / (
            f"{args.prefix}"
            ".gene_integration_long.tsv"
        )
    )

    parameters_path = (
        args.outdir
        / (
            f"{args.prefix}"
            ".gene_ranking_parameters.tsv"
        )
    )

    config_json_path = (
        args.outdir
        / (
            f"{args.prefix}"
            ".gene_ranking_parameters.json"
        )
    )

    write_gene_integration_tsv(
        short_table,
        short_path,
        args.prefix,
    )

    write_gene_integration_tsv(
        long_table,
        long_path,
        args.prefix,
    )

    pd.DataFrame(
        flatten_config(
            config
        )
    ).to_csv(
        parameters_path,
        sep="\t",
        index=False,
    )

    with config_json_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            config,
            handle,
            indent=2,
        )

    print(
        f"Wrote: {short_path}"
    )

    print(
        f"Wrote: {long_path}"
    )

    print(
        f"Wrote: {parameters_path}"
    )

    print(
        f"Wrote: {config_json_path}"
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