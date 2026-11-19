#!/bin/sh
# Check the build against a reference output file.
# usage: tests/compare.sh mirna.fa utrs.fasta reference.txt [options]
# example: tests/compare.sh dme-miR-956-3p.fa dmel-three-prime-UTR.fasta miR-956-3p-25.txt -energy -20 -T 25
set -e
mir=$1; utr=$2; ref=$3; shift 3
tmp=$(mktemp)
"$(dirname "$0")/../mirtherm" "$mir" "$utr" "$@" > "$tmp"
if cmp -s "$tmp" "$ref"; then echo "identical to $ref"; else echo "DIFFERENT from $ref"; diff "$tmp" "$ref" | head -20; exit 1; fi
rm -f "$tmp"
