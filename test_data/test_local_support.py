#!/usr/bin/env python3
import glob
import os
import random
import re
import subprocess
import sys
import tempfile


def read_fasta(path):
    sequences = []
    with open(path) as source:
        for line in source:
            if line.startswith('>'):
                sequences.append('')
            else:
                sequences[-1] += line.strip()
    return sequences


def kmers(sequence, size):
    reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return {min(sequence[i:i + size], reverse[len(sequence) - i - size:
                                              len(sequence) - i])
            for i in range(len(sequence) - size + 1)}


def run(command):
    result = subprocess.run(command, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(result.stdout)
    return result.stdout


def main():
    megahit = os.path.abspath(sys.argv[1])
    rng = random.Random(2971)
    alphabet = 'ACGT'
    genomes = [('high', ''.join(rng.choices(alphabet, k=1800)), 1000),
               ('low', ''.join(rng.choices(alphabet, k=2500)), 15)]
    with tempfile.TemporaryDirectory() as temp_dir:
        reads = os.path.join(temp_dir, 'reads.fa')
        with open(reads, 'w') as output:
            for name, sequence, depth in genomes:
                for i in range(len(sequence) * depth // 600):
                    size = rng.randrange(400, 801)
                    start = rng.randrange(len(sequence) - size + 1)
                    bases = []
                    for base in sequence[start:start + size]:
                        value = rng.random()
                        if value < .05 / 3:
                            continue
                        if value < .1 / 3:
                            bases.append(rng.choice(alphabet))
                        elif value < .05:
                            base = rng.choice(alphabet.replace(base, ''))
                        bases.append(base)
                    output.write('>{}_{}\n{}\n'.format(name, i, ''.join(bases)))

        base_dir = os.path.join(temp_dir, 'baseline')
        filtered_dir = os.path.join(temp_dir, 'filtered')
        common = [megahit, '-r', reads, '-t', '1', '-m', '800000000',
                  '--mem-flag', '2']
        run(common + ['-o', base_dir])
        run(common + ['-o', filtered_dir, '--min-local-support', '2'])
        threaded_dir = os.path.join(temp_dir, 'threaded')
        run(common + ['-t', '2', '-o', threaded_dir, '--min-local-support', '2'])

        def first_mapping_counts(directory):
            with open(os.path.join(directory, 'log')) as source:
                match = re.search(r'Lib 0: total (\d+) reads, aligned (\d+), added (\d+) reads',
                                  source.read())
            if not match:
                raise AssertionError('Local read-mapping counts were not logged')
            return match.groups()

        if first_mapping_counts(filtered_dir) != first_mapping_counts(threaded_dir):
            raise AssertionError('Local mapping counts changed with the thread count')

        baseline = read_fasta(os.path.join(base_dir, 'final.contigs.fa'))
        filtered = read_fasta(os.path.join(filtered_dir, 'final.contigs.fa'))
        if not filtered or len(filtered) >= len(baseline):
            raise AssertionError('Support filtering did not reduce redundant output')
        observed = set()
        for sequence in filtered:
            observed.update(kmers(sequence, 21))
        for name, sequence, _ in genomes:
            expected = kmers(sequence, 21)
            if len(expected & observed) < .7 * len(expected):
                raise AssertionError('Insufficient recovery of ' + name)

        rejected = glob.glob(os.path.join(filtered_dir, 'rejected_local', 'k*.fa'))
        if not rejected or not any(read_fasta(path) for path in rejected):
            raise AssertionError('Excluded local candidates were not preserved')

        for index, (extra, message) in enumerate((
                (['--min-local-support', '-1'], 'must be between'),
                (['--min-local-support', '2147483648'], 'must be between'),
                (['--min-local-support', '2', '--no-local'], 'cannot be combined'))):
            result = subprocess.run(
                common + ['-o', os.path.join(temp_dir, 'invalid' + str(index))] + extra,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0 or message not in result.stderr:
                raise AssertionError('Invalid support options were accepted')

        for extra in ([], ['--min-local-support', '2']):
            result = subprocess.run(
                common + ['-o', filtered_dir, '--force'] + extra,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0 or 'rejected_local already exists' not in result.stderr:
                raise AssertionError('A previous candidate archive could be reused')


if __name__ == '__main__':
    main()
