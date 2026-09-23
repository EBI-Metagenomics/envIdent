#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
from collections import defaultdict
import logging

import pandas as pd

# Rank columns for BOLD and MIDORI.
_TAX_RANKS = [
    "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
]

logging.basicConfig(level=logging.DEBUG)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-t", "--taxa", required=True, type=str, help="Path to taxa file"
    )
    parser.add_argument(
        "-f", "--fwd", required=True, type=str, help="Path to DADA2 forward map file"
    )
    parser.add_argument(
        "-r", "--rev", required=False, type=str, help="Path to DADA2 reverse map file"
    )
    parser.add_argument(
        "-hd", "--headers", required=False, type=str, help="Path to fastq headers"
    )
    parser.add_argument("-s", "--sample", required=True, type=str, help="Sample ID")

    args = parser.parse_args()

    taxa = args.taxa
    fwd = args.fwd
    rev = args.rev
    headers = args.headers
    sample = args.sample

    return taxa, fwd, rev, headers, sample


def order_df(taxa_df):
    if len(taxa_df.columns) == 9:
        taxa_df = taxa_df.sort_values(_TAX_RANKS, ascending=True)
    else:
        logging.error("Data frame not the right size, something wrong.")
        exit(1)

    return taxa_df


def make_tax_assignment_dict_bold(taxa_df, asv_dict):
    tax_assignment_dict = defaultdict(int)

    for i in range(len(taxa_df)):
        sorted_index = taxa_df.index[i]
        asv_num = taxa_df.iloc[i, 0]
        asv_count = asv_dict[asv_num]

        if asv_count == 0:
            continue

        d = taxa_df.loc[sorted_index, "Domain"]
        k = taxa_df.loc[sorted_index, "Kingdom"]
        p = taxa_df.loc[sorted_index, "Phylum"]
        c = taxa_df.loc[sorted_index, "Class"]
        o = taxa_df.loc[sorted_index, "Order"]
        f = taxa_df.loc[sorted_index, "Family"]
        g = taxa_df.loc[sorted_index, "Genus"]
        s = taxa_df.loc[sorted_index, "Species"]

        tax_assignment = ""

        while True:

            if d != "0":
                d = "_".join(d.split(" "))
                tax_assignment += d
            else:
                break

            if k != "0":
                k = "_".join(k.split(" "))
                tax_assignment += f"\t{k}"
            elif d != "0":
                tax_assignment += "\tk__"
            else:
                break

            if p != "0":
                p = "_".join(p.split(" "))
                tax_assignment += f"\t{p}"
            else:
                break

            if c != "0":
                c = "_".join(c.split(" "))
                tax_assignment += f"\t{c}"
            else:
                break

            if o != "0":
                o = "_".join(o.split(" "))
                tax_assignment += f"\t{o}"
            else:
                break

            if f != "0":
                f = "_".join(f.split(" "))
                tax_assignment += f"\t{f}"
            else:
                break

            if g != "0":
                g = "_".join(g.split(" "))
                tax_assignment += f"\t{g}"
            else:
                break

            if s != "0":
                s = "_".join(s.split(" "))
                tax_assignment += f"\t{s}"
            break

        if tax_assignment == "":
            continue

        tax_assignment_dict[tax_assignment] += asv_count

    return tax_assignment_dict


def generate_asv_count_dict(asv_dict):

    res_dict = defaultdict(list)

    for asv_id, count in asv_dict.items():

        if count == 0:
            continue

        res_dict["asv"].append(asv_id)
        res_dict["count"].append(count)

    res_df = pd.DataFrame(res_dict, columns=["asv", "count"])
    res_df = res_df.sort_values(by="asv", ascending=True)
    res_df = res_df.sort_values(by="count", ascending=False)

    return res_df


def main():
    taxa, fwd, rev, headers, sample = parse_args()

    fwd_fr = open(fwd, "r")
    paired_end = True

    if rev is None:
        paired_end = False
        rev_fr = [True]
    else:
        rev_fr = open(rev, "r")

    # LCA output is headerless: ASV ID followed by taxonomy, including no-hit rows.
    try:
        taxa_df = pd.read_csv(taxa, sep="\t", header=None, dtype=str)
    except pd.errors.EmptyDataError:
        taxa_df = pd.DataFrame(columns=[0, 1])
    lca = len(taxa_df.columns) == 2
    if lca:
        taxa_df.columns = ["ASV", "taxonomy"]
    else:
        # Retain support for the original headed taxonomy-table input.
        taxa_df = pd.read_csv(taxa, sep="\t", dtype=str).fillna("0")
        taxa_df = order_df(taxa_df)
    asv_list = taxa_df.ASV.to_list()

    asv_dict = defaultdict(int)

    counter = -1
    for line_fwd in fwd_fr:
        counter += 1
        line_fwd = line_fwd.strip()

        if line_fwd == "0" or f"seq_{line_fwd}" not in asv_list:
            continue

        asv_dict[f"seq_{line_fwd}"] += 1

    fwd_fr.close()
    if paired_end:
        rev_fr.close()

    if asv_dict and not lca:  # LCA/Krona aggregation is performed downstream
    
        tax_assignment_dict = make_tax_assignment_dict_bold(taxa_df,asv_dict) 

        with open(f"./{sample}_asv_krona_counts.txt", "w") as fw:
            for tax_assignment, count in tax_assignment_dict.items():
                fw.write(f"{count}\t{tax_assignment}\n")

    asv_count_df = generate_asv_count_dict(asv_dict)
    asv_count_df.to_csv(
        f"./{sample}_asv_read_counts.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
