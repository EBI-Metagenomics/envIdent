#!/usr/bin/env python

# Copyright 2025-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

"""Reformat VSEARCH taxonomy-assignment output into a fixed rank-column layout.

Reads a VSEARCH ``--userout`` TSV (query, taxonomy field, and optionally
percent identity / query coverage) and rewrites the taxonomy field as
``d__;k__;p__;c__;o__;f__;g__;s__;`` so every row has the same eight ranks
in the same order, regardless of which ranks the original hit reported.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from collections.abc import Iterable

# output rank letter -> rank word used in the input format
RANK_MAP: dict[str, str] = {
    "d": "domain",
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
WORD_TO_LETTER: dict[str, str] = {v: k for k, v in RANK_MAP.items()}
ranks: list[str] = list(RANK_MAP)  # d, k, p, c, o, f, g, s


def clean_id(val: str) -> str:
    """Strip surrounding whitespace and a trailing numeric tax ID.

    e.g. ``"Channa_gachua_64606"`` -> ``"Channa_gachua"``.
    """
    return re.sub(r"_\d+$", "", val.strip())


def parse_species(genus: str, species_raw: str, lca: bool = False) -> str:
    """Extract the species epithet from a raw ``species_...`` field value.

    Args:
        genus: The genus already parsed for this record (may be ``""``),
            used to strip a duplicated genus prefix from the species name.
        species_raw: The raw value following ``species_`` in the input,
            e.g. ``"Channa_gachua_64606"``.
        lca: If True, only the first word of the epithet is kept, and an
            epithet containing a "." (e.g. an "sp." placeholder) is
            treated as unresolved and returned as ``""``.

    Returns:
        The cleaned species epithet (empty string if ``lca`` is True and
        the epithet is unresolved).
    """
    species_name = clean_id(species_raw)

    # remove genus duplication (Channa_gachua -> gachua)
    if genus and species_name.startswith(genus + "_"):
        species_name = species_name[len(genus) + 1:]

    if lca:
        # Check the epithet after removing the matching genus prefix.
        first_word = species_name.split("_")[0]
        return "" if "." in first_word else first_word

    return species_name


def parse_tax(items: Iterable[str], lca: bool = False) -> dict[str, str]:
    """Parse ``rank_value`` entries into a rank-letter -> value mapping.

    Args:
        items: Rank entries such as ``"genus_Channa_9"``,
            ``"species_Channa_gachua_10"``. Entries whose rank word is not
            in :data:`RANK_MAP` (e.g. "clade", "subclass") are ignored.
        lca: Forwarded to :func:`parse_species` for the species entry.

    Returns:
        A dict mapping rank letters (``d``, ``k``, ``p``, ``c``, ``o``,
        ``f``, ``g``, ``s``) to their cleaned values. Ranks absent from
        ``items`` are simply absent from the returned dict.
    """
    rank_map: dict[str, str] = {}
    genus = ""

    for item in items:
        item = item.strip()
        if "_" not in item:
            continue

        rank_word, val = item.split("_", 1)
        letter = WORD_TO_LETTER.get(rank_word)
        if letter is None:
            continue  # ignores clade, subclass, etc.

        if letter == "g":
            genus = clean_id(val)
            rank_map["g"] = genus
        elif letter == "s":
            rank_map["s"] = parse_species(genus, val, lca=lca)
        else:
            rank_map[letter] = clean_id(val)

    return rank_map


def build_row(
    seq_id: str,
    taxonomy: str,
    accession: str = "",
    identity: str = "",
    coverage: str = "",
    keep_accession: bool = False,
    show_identity: bool = False,
    show_coverage: bool = False,
) -> str:
    """Assemble one tab-separated output line.

    Also used to build the header row, by passing label strings (e.g.
    ``"query"``, ``"accession"``) instead of data values. The set of flags
    (``keep_accession``, ``show_identity``, ``show_coverage``) determines
    which optional columns are included, in the fixed layout:
    ``query, [accession], taxonomy, [pid], [qcov]``.

    Returns:
        The assembled line, without a trailing newline.
    """
    cols = [seq_id]
    if keep_accession:
        cols.append(accession)
    cols.append(taxonomy)
    if show_identity:
        cols.append(identity)
    if show_coverage:
        cols.append(coverage)
    return "\t".join(cols)


def process(
    line: str,
    keep_accession: bool = False,
    show_identity: bool = False,
    show_coverage: bool = False,
    lca: bool = False,
) -> str:
    """Convert one line of VSEARCH ``--userout`` TSV to the output format.

    Args:
        line: One input line: ``query<TAB>taxonomy_field[<TAB>pid[<TAB>qcov]]``.
            ``taxonomy_field`` is either ``"*"`` (no hit, written by
            VSEARCH's ``--output_no_hits``) or
            ``">accession;rank_val;rank_val;..."``.
        keep_accession: Include the hit accession column in the output.
        show_identity: Include the percent-identity column in the output.
        show_coverage: Include the query-coverage column in the output.
        lca: Forwarded to :func:`parse_tax` / :func:`parse_species`.

    Returns:
        The formatted output line, without a trailing newline.

    Raises:
        ValueError: If the line does not have at least 2 tab-separated
            columns; if ``show_identity`` is True but there is no 3rd
            column; if ``show_coverage`` is True but there is no 4th
            column; or if the taxonomy field is not ``"*"`` and does not
            have the form ``"accession;rank_val;..."``.
    """
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 2:
        raise ValueError(
            f"expected at least 2 tab-separated columns (query, taxonomy), got {len(cols)}: {line!r}"
        )
    if show_identity and len(cols) < 3:
        raise ValueError(
            "--identity/-p was requested but the line has no 3rd column "
            f"(percent identity): {line!r}"
        )
    if show_coverage and len(cols) < 4:
        raise ValueError(
            "--coverage/-q was requested but the line has no 4th column "
            f"(query coverage): {line!r}"
        )

    seq_id = cols[0]
    tax_field = cols[1]
    identity = cols[2] if len(cols) > 2 else ""  # 3rd column: percentage identity
    coverage = cols[3] if len(cols) > 3 else ""  # 4th column: query coverage

    # VSEARCH --output_no_hits uses '*' for the missing reference target.
    # Retain the query with eight empty ranks and no invented alignment metrics.
    if tax_field == "*":
        taxonomy = ";".join(f"{r}__" for r in ranks) + ";"
        return build_row(
            seq_id, taxonomy,
            identity="NA", coverage="NA",
            keep_accession=keep_accession,
            show_identity=show_identity,
            show_coverage=show_coverage,
        )

    fields = tax_field.split(";")
    if len(fields) < 2:
        raise ValueError(
            f"malformed taxonomy field, expected 'accession;rank_val;rank_val;...', got: {tax_field!r}"
        )

    # fields[0] is the hit accession (e.g. PQ519790.1.<1.>655_root_1),
    # the rest are rank entries
    accession = fields[0].lstrip(">")
    rank_map = parse_tax(fields[1:], lca=lca)
    taxonomy = ";".join(f"{r}__{rank_map.get(r, '')}" for r in ranks) + ";"

    return build_row(
        seq_id, taxonomy,
        accession=accession, identity=identity, coverage=coverage,
        keep_accession=keep_accession,
        show_identity=show_identity,
        show_coverage=show_coverage,
    )


def main() -> None:
    """CLI entry point: parse arguments and reformat the input file.

    Writes to a temporary file as each line is processed, and only makes
    the result visible at the requested destination (via an atomic rename
    for ``-o``, or by streaming to stdout) once every line has succeeded.
    On any error, nothing is written to the destination and the process
    exits with a non-zero status.
    """
    parser = argparse.ArgumentParser(
        description="Convert taxonomy headers to d__;k__;p__;... format."
    )
    parser.add_argument(
        "-i", "--input",
        type=argparse.FileType("r"),
        required=True,
        help="Input TSV file",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output TSV file (default: stdout)",
    )
    parser.add_argument(
        "-a", "--keep-accession",
        action="store_true",
        help="Include the hit accession in the output",
    )
    parser.add_argument(
        "-p", "--identity",
        action="store_true",
        help="Include percentage identity in the output",
    )
    parser.add_argument(
        "-q", "--coverage",
        action="store_true",
        help="Include query coverage in the output",
    )
    parser.add_argument(
        "--lca", action="store_true",
        help="Omit species whose first word contains a dot after stripping the genus",
    )
    args = parser.parse_args()

    def write_rows(fout) -> None:
        """Write the header and every processed row to fout, in order."""
        header = build_row(
            "query", "taxonomy",
            accession="accession", identity="pid", coverage="qcov",
            keep_accession=args.keep_accession,
            show_identity=args.identity,
            show_coverage=args.coverage,
        )
        fout.write(header + "\n")
        with args.input as fin:
            for line_num, line in enumerate(fin, start=1):
                try:
                    out = process(
                        line,
                        keep_accession=args.keep_accession,
                        show_identity=args.identity,
                        show_coverage=args.coverage,
                        lca=args.lca,
                    )
                except ValueError as err:
                    sys.exit(f"{parser.prog}: error at input line {line_num}: {err}")
                fout.write(out + "\n")

    if args.output:
        # Write to a temp file in the same directory as the requested output,
        # and only move it into place (atomic rename) once every line has
        # processed successfully. On error, the temp file is discarded and
        # the real destination is left untouched.
        dest_dir = os.path.dirname(os.path.abspath(args.output))
        tmp = tempfile.NamedTemporaryFile(
            mode="w", dir=dest_dir, delete=False, suffix=".tmp", prefix="tax_reformat_"
        )
        tmp_path: str = tmp.name
        success = False
        try:
            with tmp:
                write_rows(tmp)
            success = True
        finally:
            if not success:
                os.remove(tmp_path)
        os.chmod(tmp_path, 0o644)
        os.replace(tmp_path, args.output)  # atomic: same filesystem as dest_dir
    else:
        # Stdout can't be "un-printed" if a later line fails, and it's
        # normally piped/redirected rather than treated as a persistent
        # deliverable the way a -o file is. So for stdout we stream directly
        # (same speed as before) and accept partial output on a failure,
        # rather than paying for an extra buffer-then-copy pass.
        write_rows(sys.stdout)


if __name__ == "__main__":
    main()
