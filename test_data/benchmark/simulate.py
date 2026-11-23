#!/usr/bin/env python3
"""Generate a reproducible mixed-community, short-read-length error simulation."""
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('prefix', help='output prefix; writes .ref.fa and .reads.fa')
parser.add_argument('--error', type=float, required=True,
                    help='per-template-base probability of an error event')
parser.add_argument('--coverage', type=int, default=1000)
parser.add_argument('--real-references', action='store_true',
                    help='read flu_source.fa, corona_source.fa and herpes_source.fa from cwd')
parser.add_argument('--flu-only', action='store_true',
                    help='use the influenza references without the larger viruses')
parser.add_argument('--minor-strain-depth', type=int, default=0,
                    help='add a 2%% divergent copy of the first influenza segment at this depth')
args = parser.parse_args()
if not 0 <= args.error <= 1 or args.coverage <= 0 or args.minor_strain_depth < 0:
    parser.error('error must be in [0,1], coverage positive, and minor depth nonnegative')
if (args.flu_only or args.minor_strain_depth) and not args.real_references:
    parser.error('--flu-only and --minor-strain-depth require --real-references')

rng = random.Random(932414)
alphabet = 'ACGT'
sizes = [2341, 2341, 2233, 1778, 1565, 1413, 1027, 890]
entries = [('segment_%d' % (i + 1), size, args.coverage)
           for i, size in enumerate(sizes)]
entries += [('coronavirus_size', 30000, args.coverage),
            ('herpesvirus_size', 152000, args.coverage),
            ('low_depth_virus', 2500, 15),
            ('bacterial_background', 30000, 25)]
genomes = [(name, ''.join(rng.choices(alphabet, k=size)), depth)
           for name, size, depth in entries]
if args.real_references:
    genomes = []
    seen = {}
    paths = ('flu_source.fa',) if args.flu_only else (
        'flu_source.fa', 'corona_source.fa', 'herpes_source.fa')
    for path in paths:
        with open(path) as source:
            for record in source.read().split('>'):
                if not record:
                    continue
                lines = record.splitlines()
                name = lines[0].split()[0]
                sequence = ''.join(lines[1:]).upper().replace('R', 'A')
                if name in seen:
                    if seen[name] != sequence:
                        raise ValueError('conflicting duplicate ' + name)
                    continue
                if set(sequence) - set(alphabet):
                    raise ValueError('unsupported ambiguous bases in ' + name)
                if len(sequence) < 800:
                    raise ValueError('reference is shorter than 800 bases: ' + name)
                seen[name] = sequence
                genomes.append((name, sequence, args.coverage))
    if not genomes:
        raise ValueError('no reference sequences found')
    if args.minor_strain_depth:
        name, sequence, _ = genomes[0]
        bases = list(sequence)
        for pos in rng.sample(range(len(bases)), max(1, round(len(bases) * 0.02))):
            bases[pos] = rng.choice(alphabet.replace(bases[pos], ''))
        genomes.append((name + '_minor', ''.join(bases), args.minor_strain_depth))
    genomes += [('low_depth_virus', ''.join(rng.choices(alphabet, k=2500)), 15),
                ('bacterial_background', ''.join(rng.choices(alphabet, k=30000)), 25)]

with open(args.prefix + '.ref.fa', 'x') as ref:
    for name, seq, depth in genomes:
        ref.write('>' + name + '\n' + seq + '\n')
with open(args.prefix + '.reads.fa', 'x') as reads:
    for name, seq, depth in genomes:
        for i in range(len(seq) * depth // 600):
            size = rng.randrange(400, 801)
            start = rng.randrange(len(seq) - size + 1)
            result = []
            for base in seq[start:start + size]:
                r = rng.random()
                if r < args.error / 3:
                    continue
                if r < args.error * 2 / 3:
                    result.append(rng.choice(alphabet))
                elif r < args.error:
                    base = rng.choice(alphabet.replace(base, ''))
                result.append(base)
            reads.write('>%s_%d\n%s\n' % (name, i, ''.join(result)))
