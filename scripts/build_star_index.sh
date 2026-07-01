#!/usr/bin/env bash

set -euo pipefail

FASTA="$1"
GTF="$2"
OUTDIR="$3"

STAR \
  --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir "$OUTDIR" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 149