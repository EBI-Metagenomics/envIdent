#!/usr/bin/env node
const fs = require('fs');

// Keep every ASV count, including ASVs with no reference hit.
const [countsFile, assignmentsFile, outputFile, countedAssignmentsFile] = process.argv.slice(2);
const assignments = new Map();
const originalTaxonomy = new Map();
for (const line of fs.readFileSync(assignmentsFile, 'utf8').split(/\r?\n/).filter(Boolean)) {
    const [id, taxonomy] = line.split('\t');
    if (assignments.has(id)) throw new Error(`Duplicate LCA assignment: ${id}`);
    originalTaxonomy.set(id, taxonomy);
    const ranks = taxonomy.split(';').filter(Boolean);
    // Stop at the deepest assigned rank, retaining gaps within the lineage.
    while (ranks.length && ranks[ranks.length - 1].endsWith('__')) ranks.pop();
    assignments.set(id, ranks.length ? ranks.join('\t') : 'Unclassified');
}
const lines = fs.readFileSync(countsFile, 'utf8').trim().split(/\r?\n/);
if (lines.shift() !== 'asv\tcount') throw new Error('Expected ASV count header: asv<TAB>count');
const totals = new Map();
const seen = new Set();
const asvCounts = new Map();
for (const line of lines.filter(Boolean)) {
    const [id, value] = line.split('\t');
    const count = Number(value);
    if (!value || !Number.isSafeInteger(count) || count < 0 || seen.has(id)) {
        throw new Error(`Invalid or duplicate ASV count: ${line}`);
    }
    seen.add(id);
    asvCounts.set(id, count);
    if (count === 0) continue;
    const taxonomy = assignments.get(id) || 'Unclassified';
    totals.set(taxonomy, (totals.get(taxonomy) || 0) + count);
}
fs.writeFileSync(outputFile, [...totals].map(([taxonomy, count]) => `${count}\t${taxonomy}\n`).join(''));

// Preserve LCA row order and taxonomy, appending the shared map-derived count.
if (countedAssignmentsFile) {
    fs.writeFileSync(countedAssignmentsFile, [...originalTaxonomy]
        .map(([id, taxonomy]) => `${id}\t${taxonomy}\t${asvCounts.get(id) || 0}\n`).join(''));
}
