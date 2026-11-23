#!/usr/bin/env python3
import os
import random
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


def canonical(sequence):
    reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return min(sequence, reverse)


def main():
    megahit = os.path.abspath(sys.argv[1])
    random_gen = random.Random(72)
    sequences = {
        size: ''.join(random_gen.choices('ACGT', k=size))
        for size in (180, 250, 400)
    }

    with tempfile.TemporaryDirectory() as temp_dir:
        reads = os.path.join(temp_dir, 'reads.fa')
        with open(reads, 'w') as output:
            for size, sequence in sequences.items():
                for i in range(50):
                    output.write('>{}_{}\n{}\n'.format(size, i, sequence))

        cases = ((200, False, (400,)), (200, True, (250, 400)),
                 (251, True, (400,)))
        for minimum, retain, expected in cases:
            out_dir = os.path.join(temp_dir, 'out{}_{}'.format(minimum, retain))
            extra = ['--retain-short-standalone'] if retain else []
            subprocess.run(
                [megahit, '-r', reads, '-o', out_dir, '-t', '2',
                 '-m', '1000000000', '--no-local',
                 '--min-contig-len', str(minimum)] + extra,
                check=True, stdout=subprocess.DEVNULL)
            actual = sorted(canonical(s) for s in
                            read_fasta(os.path.join(out_dir, 'final.contigs.fa')))
            assert actual == sorted(canonical(sequences[size]) for size in expected)

        rejected = subprocess.run(
            [megahit, '-r', reads, '-o', os.path.join(temp_dir, 'invalid'),
             '--retain-short-standalone'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert rejected.returncode != 0
        assert '--retain-short-standalone requires --no-local' in rejected.stderr


if __name__ == '__main__':
    main()
