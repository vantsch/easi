#!/usr/bin/env python3

from pathlib import Path

import pandas as pd


EVENT_FILES = {
    "alt3prime",
    "alt3and5prime",
    "alt5prime",
    "alternate_first_exon",
    "alternate_last_exon",
    "alternative_intron",
    "cassette",
    "junctions",
    "multi_exon_spanning",
    "mutually_exclusive",
    "orphan_junction",
    "tandem_cassette",
}

SKIP_FILES = {
    "heatmap",
    "other",
    "p_alt3prime",
    "p_alt5prime",
    "p_alternate_first_exon",
    "p_alternate_last_exon",
}

SHEET_NAMES = {
    "alt3prime": "Alternative 3' SS",
    "alt3and5prime": "Alternative 3' and 5' SS",
    "alt5prime": "Alternative 5' SS",
    "alternate_first_exon": "Alternative First Exon",
    "alternate_last_exon": "Alternative Last Exon",
    "alternative_intron": "Alternative Intron",
    "cassette": "Exon Skipping",
    "junctions": "Junctions",
    "multi_exon_spanning": "Multi-Exon Spanning",
    "mutually_exclusive": "Mutually Exclusive Exons",
    "orphan_junction": "Orphan Junctions",
    "tandem_cassette": "Tandem Exon Skipping",
}

SUMMARY_COLUMNS = [
    "gene_name",
    "gene_id",
    "seqid",
    "strand",
    "module_id",
    "module_event_combination",
    "num_events",
    "denovo_juncs",
]


def add_chr_prefix(seqid):
    """Return a chromosome name with a chr prefix."""

    seqid = str(seqid).strip()

    if not seqid or seqid.lower() == "nan":
        return ""

    if not seqid.startswith("chr"):
        seqid = f"chr{seqid}"

    return seqid


def format_coordinates(rep):
    """Return chromosome-qualified coordinates from the MAJIQ LSV ID."""

    lsv_id = str(rep.get("lsv_id", "")).strip()

    if not lsv_id or lsv_id.lower() == "nan":
        return ""

    coordinate = lsv_id.rsplit(":", 1)[-1].strip()

    if not coordinate:
        return ""

    seqid = add_chr_prefix(rep.get("seqid", ""))

    return (
        f"{seqid}:{coordinate}"
        if seqid
        else coordinate
    )


def format_junction_coordinates(rep):
    """Return chromosome-qualified junction coordinates."""

    coordinate = str(
        rep.get("junction_coord", "")
    ).strip()

    if not coordinate or coordinate.lower() == "nan":
        return ""

    if ":" in coordinate:
        return (
            coordinate
            if coordinate.startswith("chr")
            else f"chr{coordinate}"
        )

    seqid = add_chr_prefix(rep.get("seqid", ""))

    return (
        f"{seqid}:{coordinate}"
        if seqid
        else coordinate
    )


def choose_representative(group, probability_col):
    """Select the junction with the highest probability of changing."""

    valid = group.dropna(
        subset=[probability_col]
    )

    if valid.empty:
        return None

    return valid.loc[
        valid[probability_col].idxmax()
    ]


def describe_exon_skipping(rep, dpsi):
    """Describe cassette-exon inclusion or skipping."""

    spliced_with = str(
        rep.get("spliced_with", "")
    ).strip()

    if spliced_with == "A":
        return (
            "More exon inclusion"
            if dpsi > 0
            else "More exon skipping"
        )

    if spliced_with in {"C1", "C2"}:
        return (
            "More exon skipping"
            if dpsi > 0
            else "More exon inclusion"
        )

    return (
        "Increased exon inclusion"
        if dpsi > 0
        else "Decreased exon inclusion"
    )


def describe_change(rep, event_type, dpsi_col):
    """Translate a selected MAJIQ junction into readable wording."""

    event_type = event_type.lower()
    dpsi = float(rep[dpsi_col])

    if dpsi == 0:
        return "No directional change"

    if event_type == "cassette":
        return describe_exon_skipping(
            rep,
            dpsi,
        )

    if event_type == "tandem_cassette":
        return (
            "Increased tandem-exon junction usage"
            if dpsi > 0
            else "Decreased tandem-exon junction usage"
        )

    if event_type == "alt3prime":
        return (
            "Increased alternative 3' splice-site usage"
            if dpsi > 0
            else "Decreased alternative 3' splice-site usage"
        )

    if event_type == "alt5prime":
        return (
            "Increased alternative 5' splice-site usage"
            if dpsi > 0
            else "Decreased alternative 5' splice-site usage"
        )

    if event_type == "alt3and5prime":
        return (
            "Increased alternative splice-site usage"
            if dpsi > 0
            else "Decreased alternative splice-site usage"
        )

    if event_type == "alternative_intron":
        return (
            "Increased intron-event junction usage"
            if dpsi > 0
            else "Decreased intron-event junction usage"
        )

    if event_type == "mutually_exclusive":
        return (
            "Shift toward selected exon junction"
            if dpsi > 0
            else "Shift away from selected exon junction"
        )

    if event_type == "alternate_first_exon":
        return (
            "Increased first-exon junction usage"
            if dpsi > 0
            else "Decreased first-exon junction usage"
        )

    if event_type == "alternate_last_exon":
        return (
            "Increased last-exon junction usage"
            if dpsi > 0
            else "Decreased last-exon junction usage"
        )

    if event_type in {
        "junctions",
        "orphan_junction",
        "multi_exon_spanning",
    }:
        return (
            "Increased junction usage"
            if dpsi > 0
            else "Decreased junction usage"
        )

    return (
        "Increased event usage"
        if dpsi > 0
        else "Decreased event usage"
    )


def summarise(df, event_type):
    """Collapse each MAJIQ event ID to one representative junction."""

    psi_cols = [
        col
        for col in df.columns
        if col.endswith("_median_psi")
    ]

    dpsi_cols = [
        col
        for col in df.columns
        if col.endswith("_median_dpsi")
    ]

    probability_cols = [
        col
        for col in df.columns
        if col.endswith(
            "_probability_changing"
        )
    ]

    if (
        "event_id" not in df.columns
        or len(psi_cols) < 2
        or not dpsi_cols
        or not probability_cols
    ):
        return pd.DataFrame()

    dpsi_col = dpsi_cols[0]
    probability_col = probability_cols[0]

    rows = []

    for event_id, group in df.groupby(
        "event_id",
        dropna=False,
    ):
        rep = choose_representative(
            group,
            probability_col,
        )

        if rep is None:
            continue

        dpsi = float(rep[dpsi_col])
        probability = float(
            rep[probability_col]
        )

        row = {
            "Gene": rep.get(
                "gene_name",
                "",
            ),
            "Gene ID": rep.get(
                "gene_id",
                "",
            ),
            "Chr": add_chr_prefix(
                rep.get("seqid", "")
            ),
            "Strand": rep.get(
                "strand",
                "",
            ),
            "Coordinates": format_coordinates(
                rep
            ),
            "Junction coordinates":
                format_junction_coordinates(rep),
            "Event ID": event_id,
            "Change": describe_change(
                rep,
                event_type,
                dpsi_col,
            ),
            psi_cols[0]: round(
                float(rep[psi_cols[0]]),
                3,
            ),
            psi_cols[1]: round(
                float(rep[psi_cols[1]]),
                3,
            ),
            "dPSI": round(
                dpsi,
                3,
            ),
            "|ΔPSI|": round(
                abs(dpsi),
                3,
            ),
            "Significant": (
                "Yes"
                if probability >= 0.95
                else "No"
            ),
        }

        rows.append(row)

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    return (
        out.sort_values(
            [
                "Significant",
                "|ΔPSI|",
                "Gene",
            ],
            ascending=[
                False,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )


def prepare_summary(df):
    """Select and clarify the main MAJIQ summary columns."""

    columns = [
        col
        for col in SUMMARY_COLUMNS
        if col in df.columns
    ]

    if not columns:
        return pd.DataFrame()

    out = df[columns].copy()

    if "module_event_combination" in out.columns:
        out["module_event_combination"] = (
            out["module_event_combination"]
            .astype(str)
            .str.replace(
                "tandem_cassette",
                "tandem exon skipping",
                regex=False,
            )
            .str.replace(
                "cassette",
                "exon skipping",
                regex=False,
            )
        )

    return out


def collapse_duplicate_events(events):
    """
    Collapse events that occur in more than one MAJIQ event table.

    The generic junctions table can overlap with classified event
    tables such as cassette or alt3prime. Each Event ID is counted
    only once in the gene-level summary.
    """

    rows = []

    for (
        gene_id,
        gene_name,
        event_id,
    ), group in events.groupby(
        [
            "Gene ID",
            "Gene",
            "Event ID",
        ],
        dropna=False,
    ):
        effects = pd.to_numeric(
            group["dPSI"],
            errors="coerce",
        )

        abs_effects = effects.abs()

        if abs_effects.notna().any():
            representative_index = (
                abs_effects.idxmax()
            )
            representative_effect = effects.loc[
                representative_index
            ]
            max_abs_effect = abs_effects.max()
        else:
            representative_effect = 0.0
            max_abs_effect = 0.0

        event_types = set(
            group["event_type"]
            .dropna()
            .astype(str)
        )

        # Remove the generic label when a more informative
        # event classification is available.
        if (
            "junctions" in event_types
            and len(event_types) > 1
        ):
            event_types.remove("junctions")

        rows.append(
            {
                "Gene ID": gene_id,
                "Gene": gene_name,
                "Event ID": event_id,
                "significant": bool(
                    group["is_significant"].any()
                ),
                "dPSI": representative_effect,
                "max_abs_effect": max_abs_effect,
                "event_types": "|".join(
                    sorted(event_types)
                ),
            }
        )

    return pd.DataFrame(rows)


def determine_effect_direction(group):
    """
    Determine the direction of MAJIQ effects for one gene.

    Significant events are used when present. Otherwise, all detected
    events are used.
    """

    significant = group[
        group["significant"]
    ]

    direction_group = (
        significant
        if not significant.empty
        else group
    )

    effects = pd.to_numeric(
        direction_group["dPSI"],
        errors="coerce",
    ).dropna()

    effects = effects[
        effects != 0
    ]

    if effects.empty:
        return "none"

    has_positive = bool(
        (effects > 0).any()
    )
    has_negative = bool(
        (effects < 0).any()
    )

    if has_positive and has_negative:
        return "mixed"

    if has_positive:
        return "increased"

    return "decreased"


def prepare_gene_summary(event_tables):
    """Create one standardized gene-level MAJIQ summary."""

    columns = [
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

    if not event_tables:
        return pd.DataFrame(
            columns=columns
        )

    events = pd.concat(
        event_tables,
        ignore_index=True,
    )

    if events.empty:
        return pd.DataFrame(
            columns=columns
        )

    events["is_significant"] = (
        events["Significant"]
        .astype(str)
        .str.lower()
        .eq("yes")
    )

    collapsed = collapse_duplicate_events(
        events
    )

    if collapsed.empty:
        return pd.DataFrame(
            columns=columns
        )

    rows = []

    for (
        gene_id,
        gene_name,
    ), group in collapsed.groupby(
        [
            "Gene ID",
            "Gene",
        ],
        dropna=False,
    ):
        event_types = set()

        for value in (
            group["event_types"]
            .dropna()
            .astype(str)
        ):
            event_types.update(
                item
                for item in value.split("|")
                if item
            )

        rows.append(
            {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "method": "MAJIQ",
                "detected": True,
                "significant": bool(
                    group["significant"].any()
                ),
                "n_events": int(
                    len(group)
                ),
                "n_significant_events": int(
                    group["significant"].sum()
                ),
                "max_abs_effect": round(
                    pd.to_numeric(
                        group["max_abs_effect"],
                        errors="coerce",
                    ).max(),
                    3,
                ),
                "effect_direction":
                    determine_effect_direction(
                        group
                    ),
                "event_types": "|".join(
                    sorted(event_types)
                ),
            }
        )

    result = pd.DataFrame(
        rows,
        columns=columns,
    )

    return (
        result.sort_values(
            [
                "significant",
                "n_significant_events",
                "max_abs_effect",
                "gene_name",
            ],
            ascending=[
                False,
                False,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )


def main():
    output_file = Path(
        "MAJIQ_summary.xlsx"
    )

    gene_summary_file = Path(
        "MAJIQ_gene_summary.tsv"
    )

    writer = pd.ExcelWriter(
        output_file,
        engine="openpyxl",
    )

    sheets_written = 0
    event_tables = []

    for file in sorted(
        Path(".").glob("*.tsv")
    ):
        print(
            f"Processing {file.name}"
        )

        df = pd.read_csv(
            file,
            sep="\t",
            comment="#",
        )

        if df.empty:
            print(
                f"Skipping {file.name}: "
                "empty table"
            )
            continue

        name = file.stem.lower()

        if name == "summary":
            summary = prepare_summary(
                df
            )

            if summary.empty:
                print(
                    f"Skipping {file.name}: "
                    "no expected summary columns"
                )
                continue

            summary.to_excel(
                writer,
                sheet_name="Summary",
                index=False,
            )

            sheets_written += 1
            continue

        if name in SKIP_FILES:
            print(
                f"Skipping {file.name}: "
                "excluded output type"
            )
            continue

        if name not in EVENT_FILES:
            print(
                f"Skipping {file.name}: "
                "unsupported output type"
            )
            continue

        if "event_id" not in df.columns:
            print(
                f"Skipping {file.name}: "
                "no event_id column"
            )
            continue

        out = summarise(
            df,
            name,
        )

        if out.empty:
            print(
                f"Skipping {file.name}: "
                "no reportable events"
            )
            continue

        out["event_type"] = name

        event_tables.append(
            out.copy()
        )

        sheet_name = SHEET_NAMES.get(
            name,
            file.stem
            .replace("_", " ")
            .title(),
        )[:31]

        out.drop(
            columns=["event_type"],
        ).to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )

        sheets_written += 1

    if sheets_written == 0:
        writer.close()

        output_file.unlink(
            missing_ok=True
        )

        raise RuntimeError(
            "No worksheets were generated "
            "from the MAJIQ TSV files."
        )

    writer.close()

    gene_summary = prepare_gene_summary(
        event_tables
    )

    gene_summary.to_csv(
        gene_summary_file,
        sep="\t",
        index=False,
    )

    print(
        f"Finished: {output_file}"
    )

    print(
        f"Finished: {gene_summary_file}"
    )


if __name__ == "__main__":
    main()