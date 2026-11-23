#!/usr/bin/env python3
import os
import subprocess
import sys
import tempfile


def main():
    evaluator = os.path.join(os.path.dirname(__file__), 'benchmark', 'evaluate.py')
    with tempfile.TemporaryDirectory() as directory:
        reference = os.path.join(directory, 'reference.fa')
        contigs = os.path.join(directory, 'contigs.fa')
        aligner = os.path.join(directory, 'aligner')
        with open(reference, 'w') as output:
            output.write('>r1\n' + 'A' * 2000 + '\n>r2\n' + 'A' * 300 + '\n')
        with open(contigs, 'w') as output:
            output.write('>c1\n' + 'A' * 500 + '\n>c2\n' + 'A' * 400 + '\n')
        lines = [
            'c1\t500\t0\t200\t+\tr1\t2000\t0\t200\t190\t200\t60\ttp:A:P',
            'c1\t500\t200\t400\t+\tr1\t2000\t200\t400\t190\t200\t60\ttp:A:P',
            'c2\t400\t0\t100\t+\tr1\t2000\t100\t200\t95\t100\t60\ttp:A:P',
        ]
        with open(aligner, 'w') as output:
            output.write('#!' + sys.executable + '\nprint(' + repr('\n'.join(lines)) + ')\n')
        os.chmod(aligner, 0o755)
        result = subprocess.run(
            [sys.executable, evaluator, aligner, reference, contigs],
            check=True, stdout=subprocess.PIPE, text=True)
        rows = {line.split()[0]: line.split()[1:]
                for line in result.stdout.splitlines()[2:]}
        metrics = dict(zip(rows['r1'][::2], rows['r1'][1::2]))
        assert metrics['covered'] == '400'
        assert metrics['alignments'] == '3'
        assert metrics['redundant_ref_bases'] == '100'
        assert metrics['contigs'] == '2'
        assert metrics['contig_bases'] == '900'
        metrics = dict(zip(rows['r2'][::2], rows['r2'][1::2]))
        assert metrics['contigs'] == '0' and metrics['contig_bases'] == '0'


if __name__ == '__main__':
    main()
