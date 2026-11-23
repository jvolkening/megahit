#!/usr/bin/env python3
import collections
import os
import subprocess
import sys
import tempfile


def main():
    simulator = os.path.join(os.path.dirname(__file__), 'benchmark', 'simulate.py')
    with tempfile.TemporaryDirectory() as directory:
        with open(os.path.join(directory, 'flu_source.fa'), 'w') as source:
            source.write('>segment\n' + 'ACGT' * 250 + '\n')

        for mode, extra in (('synthetic', []),
                            ('reference', ['--real-references', '--flu-only'])):
            for name, options, depths in (
                    ('default', [], (15, 25)),
                    ('custom', ['--virus-control-depth', '10',
                                '--background-depth', '10'], (10, 10))):
                prefix = os.path.join(directory, mode + '-' + name)
                subprocess.run([sys.executable, simulator, prefix, '--coverage', '1',
                                '--error', '0.02'] + extra + options,
                               cwd=directory, check=True)
                counts = collections.Counter()
                with open(prefix + '.reads.fa') as reads:
                    for line in reads:
                        if line.startswith('>low_depth_virus_'):
                            counts['virus'] += 1
                        elif line.startswith('>bacterial_background_'):
                            counts['background'] += 1
                assert counts['virus'] == 2500 * depths[0] // 600
                assert counts['background'] == 30000 * depths[1] // 600

        for option in ('--virus-control-depth', '--background-depth'):
            result = subprocess.run(
                [sys.executable, simulator, os.path.join(directory, 'invalid'),
                 '--error', '0.02', option, '0'],
                cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True)
            assert result.returncode != 0 and 'control depths positive' in result.stderr


if __name__ == '__main__':
    main()
