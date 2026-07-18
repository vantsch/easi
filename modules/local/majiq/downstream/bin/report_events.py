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
    """
    Return chromosome-qualified coordinates from the MAJIQ LSV ID.

    Example:
        lsv_id = ENSG00000182944:s:29273738-29273864
        seqid  = 22

    Returns:
        chr22:29273738-29273864
    """

    lsv_id = str(rep.get("lsv_id", "")).strip()

    if not lsv_id or lsv_id.lower() == "nan":
        return ""

    coordinate = lsv_id.rsplit(":", 1)[-1].strip()

    if not coordinate:
        return ""

    seqid = add_chr_prefix(rep.get("seqid", ""))

    if not seqid:
        return coordinate

    return f"{seqid}:{coordinate}"


def format_junction_coordinates(rep):
    """Return chromosome-qualified junction coordinates."""

    coordinate = str(rep.get("junction_coord", "")).strip()

    if not coordinate or coordinate.lower() == "nan":
        return ""

    if ":" in coordinate:
        return coordinate if coordinate.startswith("chr") else f"chr{coordinate}"

    seqid = add_chr_prefix(rep.get("seqid", ""))

    return f"{seqid}:{coordinate}" if seqid else coordinate


def choose_representative(group, probability_col):
    """Select the junction with the highest probability of changing."""

    valid = group.dropna(subset=[probability_col])

    if valid.empty:
        return None

    return valid.loc[valid[probability_col].idxmax()]


def describe_exon_skipping(rep, dpsi):
    """Describe cassette-exon inclusion or skipping."""

    spliced_with = str(rep.get("spliced_with", "")).strip()

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
        return describe_exon_skipping(rep, dpsi)

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
    """Collapse each MAJIQ event_id to one representative junction."""

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
        if col.endswith("_probability_changing")
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

        row = {
            "Gene": rep.get("gene_name", ""),
            "Gene ID": rep.get("gene_id", ""),
            "Chr": add_chr_prefix(rep.get("seqid", "")),
            "Strand": rep.get("strand", ""),
            "Coordinates": format_coordinates(rep),
            "Junction coordinates": format_junction_coordinates(rep),
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
            "|ΔPSI|": round(
                abs(float(rep[dpsi_col])),
                3,
            ),
            "Significant": (
                "Yes"
                if float(rep[probability_col]) >= 0.95
                else "No"
            ),
        }

        rows.append(row)

    out = pd.DataFrame(rows)

    if out.empty:
        return out

    return out.sort_values(
        ["Significant", "|ΔPSI|", "Gene"],
        ascending=[False, False, True],
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


def main():
    output_file = Path("MAJIQ_summary.xlsx")

    writer = pd.ExcelWriter(
        output_file,
        engine="openpyxl",
    )

    sheets_written = 0

    for file in sorted(Path(".").glob("*.tsv")):
        print(f"Processing {file.name}")

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
            summary = prepare_summary(df)

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

        sheet_name = SHEET_NAMES.get(
            name,
            file.stem.replace("_", " ").title(),
        )[:31]

        out.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False,
        )

        sheets_written += 1

    if sheets_written == 0:
        writer.close()
        output_file.unlink(missing_ok=True)

        raise RuntimeError(
            "No worksheets were generated "
            "from the MAJIQ TSV files."
        )

    writer.close()

    print(f"Finished: {output_file}")


if __name__ == "__main__":
    main()