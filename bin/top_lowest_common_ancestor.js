#!/usr/bin/env node

/*
 * Copyright 2025-2026 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 */

const fs = require("fs");
const readline = require("readline");

// Rank order from finest to coarsest, with depth = number of lineage
// components to keep (sk, k, p, c, o, f, g, s -> depths 1..8) and
// default PID / qcov thresholds. Override via CLI args.
const RANKS = [
  { name: "species",      depth: 8, pidThreshold: 96, qcovThreshold: 0 },
  { name: "genus",        depth: 7, pidThreshold: 88, qcovThreshold: 0 },
  { name: "family",       depth: 6, pidThreshold: 84, qcovThreshold: 0 },
  { name: "order",        depth: 5, pidThreshold: 81, qcovThreshold: 0 },
  { name: "class",        depth: 4, pidThreshold: 78, qcovThreshold: 0 },
  { name: "phylum",       depth: 3, pidThreshold: 77, qcovThreshold: 0 },
  { name: "kingdom",      depth: 2, pidThreshold: 77,  qcovThreshold: 0 },
  { name: "superkingdom", depth: 1, pidThreshold: 77,  qcovThreshold: 0  },
];

// Rank prefixes in order from depth 1 (superkingdom) to depth 8 (species).
// Used to pad truncated lineages with empty placeholders for dropped ranks.
const PREFIXES = ["sk__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"];

// Parse "species=95,genus=92,family=85" style overrides.
function parseThresholdArg(str) {
  const overrides = {};
  if (!str) return overrides;
  for (const pair of str.split(",")) {
    const [name, val] = pair.split("=").map(s => s.trim());
    if (name && val !== undefined && !Number.isNaN(Number(val))) {
      overrides[name.toLowerCase()] = Number(val);
    }
  }
  return overrides;
}

function buildRanks(pidOverrideStr, qcovOverrideStr) {
  const pidOverrides = parseThresholdArg(pidOverrideStr);
  const qcovOverrides = parseThresholdArg(qcovOverrideStr);
  return RANKS.map(r => ({
    ...r,
    pidThreshold: pidOverrides[r.name] !== undefined ? pidOverrides[r.name] : r.pidThreshold,
    qcovThreshold: qcovOverrides[r.name] !== undefined ? qcovOverrides[r.name] : r.qcovThreshold,
  }));
}

function parseLineage(str) {
  const entries = str.split(";").map(x => x.trim()).filter(Boolean);
  return PREFIXES.map(prefix => entries.find(entry => entry.startsWith(prefix)) || prefix);
}

function lca(paths) {
  if (paths.length === 0) return [];

  const base = paths[0];
  let len = base.length;

  for (let i = 1; i < paths.length; i++) {
    const cur = paths[i];
    let j = 0;

    while (j < len && j < cur.length && base[j] === cur[j]) {
      j++;
    }

    len = j;
    if (len === 0) break;
  }

  return base.slice(0, len);
}

// Pad a truncated lineage array out to the full 8 ranks, filling any
// dropped ranks with their empty prefix (e.g. "f__", "g__", "s__").
function padLineage(lineage) {
  const padded = lineage.slice();
  for (let i = padded.length; i < PREFIXES.length; i++) {
    padded.push(PREFIXES[i]);
  }
  return padded;
}

// Walk ranks finest -> coarsest. At the first rank where some hits clear
// BOTH the PID and qcov thresholds, narrow to just the hits with the
// highest PID among those passing, compute their LCA, and truncate it
// to that rank's depth, then stop.
function rankedLca(hits, ranks) {
  for (const rank of ranks) {
    const passing = hits.filter(
      h => h.pid >= rank.pidThreshold && h.qcov >= rank.qcovThreshold
    );
    if (passing.length === 0) continue;

    const maxPid = Math.max(...passing.map(h => h.pid));
    const topHits = passing.filter(h => h.pid === maxPid);
    const maxQcov = Math.max(...topHits.map(h => h.qcov));

    const common = lca(topHits.map(h => h.lineage));
    return { lineage: common.slice(0, rank.depth), pid: maxPid, qcov: maxQcov };
  }
  return { lineage: [], pid: null, qcov: null }; // nothing passed even the coarsest thresholds
}

async function main(file, pidThresholdArg, qcovThresholdArg) {
  const ranks = buildRanks(pidThresholdArg, qcovThresholdArg);
  const groups = new Map();

  const rl = readline.createInterface({
    input: fs.createReadStream(file),
    crlfDelay: Infinity,
  });

  for await (const line of rl) {
    if (!line.trim()) continue;

    const [id, lineage, pidStr, qcovStr] = line.split(/\t+/);
    if (!id || !lineage || (id === "query" && lineage === "taxonomy")) continue;

    const pid = parseFloat(pidStr);
    const qcov = parseFloat(qcovStr);
    const parsed = parseLineage(lineage);

    if (!groups.has(id)) groups.set(id, []);
    groups.get(id).push({ lineage: parsed, pid, qcov });
  }

  for (const [id, hits] of groups) {
    const { lineage: common, pid, qcov } = rankedLca(hits, ranks);
    const padded = padLineage(common);
    console.log(`${id}\t${padded.join(";")};`);
  }
}

// Usage:
//   node lca.js input.tsv "species=95,genus=92,family=85" "species=90,genus=85,family=80"
//                          ^-- PID overrides                ^-- qcov overrides
main(process.argv[2], process.argv[3], process.argv[4]);
