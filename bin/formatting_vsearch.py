#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

import sys
import re
import argparse

# output rank letter -> rank word used in the input format
RANK_MAP = {
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
WORD_TO_LETTER = {v: k for k, v in RANK_MAP.items()}
ranks = list(RANK_MAP)  # k, p, c, o, f, g, s


def clean_id(val):
    # remove trailing numeric tax IDs like _64606
    return re.sub(r"_\d+$", "", val.strip())


def parse_species(genus, species_raw):
    species_name = clean_id(species_raw)

    # remove genus duplication (Channa_gachua -> gachua)
    if genus and species_name.startswith(genus + "_"):
        species_name = species_name[len(genus) + 1:]

    # fallback: keep only epithet
    if "_" in species_name:
        species_name = species_name.split("_")[0]

    return species_name


def parse_tax(items):
    rank_map = {}
    genus = ""

    for item in items:
        item = item.strip()
        if "_" not in item:
            continue

        rank_word, val = item.split("_", 1)
        letter = WORD_TO_LETTER.get(rank_word)
        if letter is None:
            continue  # ignores domain, clade, subclass, etc.

        if letter == "g":
            genus = clean_id(val)
            rank_map["g"] = genus
        elif letter == "s":
            rank_map["s"] = parse_species(genus, val)
        else:
            rank_map[letter] = clean_id(val)

    return rank_map


def make_header(keep_accession=False, show_identity=False, show_coverage=False):
    header = ["query"]
    if keep_accession:
        header.append("accession")
    header.append("taxonomy")
    if show_identity:
        header.append("pid")
    if show_coverage:
        header.append("qcov")
    return "\t".join(header)


def process(line, keep_accession=False, show_identity=False, show_coverage=False):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 2:
        return None

    seq_id = cols[0]
    tax_field = cols[1]
    identity = cols[2] if len(cols) > 2 else ""  # 3rd column: percentage identity
    coverage = cols[3] if len(cols) > 3 else ""  # 4th column: query coverage

    # VSEARCH --output_no_hits uses '*' for the missing reference target.
    # Retain the query with eight empty ranks and no invented alignment metrics.
    if tax_field == "*":
        out_cols = [seq_id]
        if keep_accession:
            out_cols.append("")
        out_cols.append("k__;p__;c__;o__;f__;g__;s__;")
        if show_identity:
            out_cols.append("NA")
        if show_coverage:
            out_cols.append("NA")
        return "\t".join(out_cols)

    fields = tax_field.split(";")
    if len(fields) < 2:
        return None

    # fields[0] is the hit accession (e.g. PQ519790.1.<1.>655_root_1),
    # the rest are rank entries
    accession = fields[0].lstrip(">")
    rank_map = parse_tax(fields[1:])

    taxonomy = ";".join(f"{r}__{rank_map.get(r, '')}" for r in ranks) + ";"
    
    out_cols = [seq_id]
    if keep_accession:
        out_cols.append(accession)
    out_cols.append(taxonomy)
    if show_identity:
        out_cols.append(identity)
    if show_coverage:
        out_cols.append(coverage)
    return "\t".join(out_cols)


def main():
    parser = argparse.ArgumentParser(
        description="Convert taxonomy headers to k__;p__;... format."
    )
    parser.add_argument(
        "-i", "--input",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input TSV file (default: stdin)",
    )
    parser.add_argument(
        "-o", "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
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
    args = parser.parse_args()

    with args.input as fin, args.output as fout:
        fout.write(
            make_header(args.keep_accession, args.identity, args.coverage) + "\n"
        )
        for line in fin:
            out = process(
                line,
                keep_accession=args.keep_accession,
                show_identity=args.identity,
                show_coverage=args.coverage,
            )
            if out:
                fout.write(out + "\n")


if __name__ == "__main__":
    main()
