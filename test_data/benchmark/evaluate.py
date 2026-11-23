#!/usr/bin/env python3
"""Summarize primary minimap2 alignments; these metrics are not a full assembly evaluation."""
import collections
import subprocess
import sys

aligner, reference, contigs = sys.argv[1:]


def fasta(path):
    with open(path) as source:
        name, seq = None, []
        for line in source:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, ''.join(seq)


def coverage(intervals):
    rightmost = total = 0
    for left, right in sorted(intervals):
        total += max(0, right - max(left, rightmost))
        rightmost = max(rightmost, right)
    return total


def lengths(path):
    records = {}
    for name, sequence in fasta(path):
        if name in records:
            raise ValueError('duplicate FASTA identifier in {}: {}'.format(path, name))
        records[name] = len(sequence)
    return records


refs = lengths(reference)
ctgs = lengths(contigs)
cmd = [aligner, '-x', 'map-ont', '-k', '11', '-w', '5', '-s', '20',
       '-c', reference, contigs]
alignments = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
targets = collections.defaultdict(list)
target_contigs = collections.defaultdict(set)
queries = collections.defaultdict(list)
query_targets = collections.defaultdict(set)
matches = blocks = 0
for line in alignments.stdout.splitlines():
    fields = line.split('\t')
    if ('tp:A:P' not in fields[12:] or int(fields[10]) < 100 or
            int(fields[9]) / int(fields[10]) < 0.8):
        continue
    targets[fields[5]].append((int(fields[7]), int(fields[8])))
    target_contigs[fields[5]].add(fields[0])
    queries[fields[0]].append((int(fields[2]), int(fields[3])))
    query_targets[fields[0]].add(fields[5])
    matches += int(fields[9])
    blocks += int(fields[10])
print('contigs', len(ctgs), 'bases', sum(ctgs.values()),
      'unaligned', sum(name not in queries for name in ctgs),
      'cross_reference', sum(len(t) > 1 for t in query_targets.values()))
print('aligned_query_bases', sum(coverage(i) for i in queries.values()),
      'alignment_identity', matches / blocks if blocks else 0)
for name, length in refs.items():
    hits = targets[name]
    covered = coverage(hits)
    redundant = sum(right - left for left, right in hits) - covered
    print(name, 'reference', length, 'covered', covered,
          'alignments', len(hits), 'redundant_ref_bases', redundant,
          'contigs', len(target_contigs[name]),
          'contig_bases', sum(ctgs[contig] for contig in target_contigs[name]))
