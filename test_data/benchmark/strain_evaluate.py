#!/usr/bin/env python3
"""Look for simulated minor-strain SNPs co-occurring in aligned contigs."""
import argparse
import bisect
import re
import subprocess


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


def unique_fasta(path):
    records = {}
    for name, sequence in fasta(path):
        if name in records:
            raise ValueError('duplicate FASTA identifier in {}: {}'.format(path, name))
        records[name] = sequence
    return records


parser = argparse.ArgumentParser()
parser.add_argument('minimap2')
parser.add_argument('reference')
parser.add_argument('contigs')
parser.add_argument('major')
parser.add_argument('minor')
args = parser.parse_args()

refs = unique_fasta(args.reference)
queries = unique_fasta(args.contigs)
major, minor = refs[args.major], refs[args.minor]
if len(major) != len(minor) or not major:
    raise ValueError('major and minor references must be nonempty and equal length')
variants = [i for i in range(len(major)) if major[i] != minor[i]]
cmd = [args.minimap2, '-x', 'map-ont', '-k', '11', '-w', '5', '-s', '20',
       '-c', '--secondary=no', args.reference, args.contigs]
alignments = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
found = []
for line in alignments.stdout.splitlines():
    fields = line.split('\t')
    if fields[5] not in (args.major, args.minor) or 'tp:A:P' not in fields[12:]:
        continue
    if int(fields[10]) < 100 or int(fields[9]) / int(fields[10]) < .8:
        continue
    seq = queries[fields[0]]
    query = int(fields[2])
    if fields[4] == '-':
        seq = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        query = int(fields[1]) - int(fields[3])
    pos = int(fields[7])
    cigar = next(field[5:] for field in fields[12:] if field.startswith('cg:Z:'))
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if ''.join(size + op for size, op in operations) != cigar:
        raise ValueError('unsupported CIGAR: ' + cigar)
    alt = major_count = 0
    seen = []
    for num, op in operations:
        size = int(num)
        if op in 'M=X':
            start = bisect.bisect_left(variants, pos)
            end = bisect.bisect_left(variants, pos + size)
            for index in variants[start:end]:
                base = seq[query + index - pos]
                if base == minor[index]:
                    alt += 1
                    seen.append(index)
                elif base == major[index]:
                    major_count += 1
            query += size
            pos += size
        elif op in 'IS':
            query += size
        elif op in 'DN':
            pos += size
    if alt >= 4 and alt > major_count:
        found.append((fields[0], alt, major_count, seen))

print('variant_sites', len(variants), 'qualifying_contigs', len(found),
      'variant_sites_in_qualifying_contigs',
      len({pos for item in found for pos in item[-1]}))
for name, alt, reference_count, _ in found:
    print(name, 'minor_alleles', alt, 'major_alleles', reference_count)
