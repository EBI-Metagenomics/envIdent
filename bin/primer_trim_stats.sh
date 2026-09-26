#!/usr/bin/env bash

# Copyright 2025-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

# Summarise cutadapt --info-file output into a table of
# pre-primer / primer / post-primer lengths and read counts.
#
# Usage: primer_trim_stats.sh <input.tsv> <output.tsv>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") <input.tsv> <output.tsv>" >&2
    exit 1
fi

input="$1"
output="$2"

printf 'Pre-primer length\tPrimer length\tPost-primer length\tRead count\n' > "${output}"

awk -F'\t' '{ print length($5), length($6), length($7) }' "${input}" \
    | sort \
    | uniq -c \
    | awk 'OFS="\t" { print $2, $3, $4, $1 }' \
    | sort -nr -k4,4 >> "${output}"
