#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


OUTPUT_COLUMNS = [
    "seq_id",
    "alignment_start",
    "alignment_end",
    "envelope_start",
    "envelope_end",
    "hmm_acc",
    "hmm_name",
    "type",
    "hmm_start",
    "hmm_end",
    "hmm_length",
    "bit_score",
    "E_value",
    "significance",
    "clan",
]


def strip_version(accession: str) -> str:
    """Convert PF00079.26 to PF00079."""
    return accession.split(".", maxsplit=1)[0]


def read_clans(path: Path) -> dict[str, str]:
    """
    Read Pfam-A.clans.tsv.

    Expected initial columns:
        Pfam accession, clan accession, ...

    Families without a clan are assigned No_clan.
    """
    clans: dict[str, str] = {}

    with path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for line_number, fields in enumerate(reader, start=1):
            if not fields or fields[0].startswith("#"):
                continue

            if len(fields) < 2:
                raise ValueError(
                    f"{path}: line {line_number} contains fewer than two columns"
                )

            pfam_accession = strip_version(fields[0].strip())
            clan_accession = fields[1].strip() or "No_clan"

            clans[pfam_accession] = clan_accession

    if not clans:
        raise ValueError(f"No Pfam clan mappings were found in {path}")

    return clans


def parse_domtblout(path: Path, clans: dict[str, str]) -> list[list[str]]:
    """
    Convert HMMER hmmscan --domtblout records to pfam_scan.pl-style records.

    HMMER domtblout fields:
      0  target name
      1  target accession
      2  target length
      3  query name
      11 conditional E-value
      12 independent E-value
      13 domain score
      15 HMM start
      16 HMM end
      17 alignment start
      18 alignment end
      19 envelope start
      20 envelope end
    """
    results: list[list[str]] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split(maxsplit=22)

            if len(fields) < 22:
                raise ValueError(
                    f"{path}: line {line_number} is not a valid "
                    f"HMMER domtblout record ({len(fields)} fields found)"
                )

            hmm_name = fields[0]
            hmm_accession = fields[1]
            hmm_length = fields[2]
            seq_id = fields[3]

            independent_evalue = fields[12]
            domain_score = fields[13]

            hmm_start = fields[15]
            hmm_end = fields[16]
            alignment_start = fields[17]
            alignment_end = fields[18]
            envelope_start = fields[19]
            envelope_end = fields[20]

            pfam_accession = strip_version(hmm_accession)
            clan = clans.get(pfam_accession, "No_clan")

            results.append(
                [
                    seq_id,
                    alignment_start,
                    alignment_end,
                    envelope_start,
                    envelope_end,
                    hmm_accession,
                    hmm_name,
                    "Domain",
                    hmm_start,
                    hmm_end,
                    hmm_length,
                    domain_score,
                    independent_evalue,
                    "1",
                    clan,
                ]
            )

    return results


def write_results(path: Path, rows: list[list[str]]) -> None:
    """Write tab-separated pfam_scan.pl-compatible output."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerows(rows)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert HMMER hmmscan domtblout output into the Pfam format "
            "accepted by IsoformSwitchAnalyzeR analyzePFAM()."
        )
    )

    parser.add_argument(
        "--domtblout",
        required=True,
        type=Path,
        help="HMMER hmmscan --domtblout file",
    )
    parser.add_argument(
        "--clans",
        required=True,
        type=Path,
        help="Decompressed Pfam-A.clans.tsv file",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV file",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    for input_path in (args.domtblout, args.clans):
        if not input_path.is_file():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")

    clans = read_clans(args.clans)
    rows = parse_domtblout(args.domtblout, clans)
    write_results(args.output, rows)

    print(
        f"Wrote {len(rows)} Pfam domain records to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)