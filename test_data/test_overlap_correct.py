#!/usr/bin/env python3
import bz2
import gzip
import os
import random
import subprocess
import sys
import tempfile


def executable(path, text):
    with open(path, 'w') as output:
        output.write('#!' + sys.executable + '\n' + text)
    os.chmod(path, 0o755)


def run(command, env, success=True, message=None):
    result = subprocess.run(command, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    if (result.returncode == 0) != success or (message and message not in result.stdout):
        raise AssertionError(result.stdout)
    return result.stdout


def main():
    megahit = os.path.abspath(sys.argv[1])
    rng = random.Random(738)
    sequence = ''.join(rng.choices('ACGT', k=500))
    with tempfile.TemporaryDirectory() as directory:
        tools = os.path.join(directory, 'tools')
        os.mkdir(tools)
        executable(os.path.join(tools, 'minimap2'),
                   'import os, sys\n'
                   'assert sys.argv[1:8] == ["-x", "ava-ont", "-k11", "-w5", "-m40", "-s40", "-c"]\n'
                   'if not os.environ.get("MOCK_NO_OVERLAPS"):\n'
                   '    print("mock overlaps")\n')
        executable(os.path.join(tools, 'racon'),
                   'import os, sys\n'
                   'assert sys.argv[1:3] == ["-f", "-u"]\n'
                   'if os.environ.get("MOCK_TRUNCATE"):\n'
                   '    print(">mh0r\\nACGT")\n'
                   '    sys.exit(0)\n'
                   'with open(sys.argv[-1]) as source:\n'
                   '    for line in source:\n'
                   '        print(line.strip() + ("r" if line.startswith(">") else ""))\n')
        env = os.environ.copy()
        env['PATH'] = tools + os.pathsep + env.get('PATH', '')
        reads = os.path.join(directory, 'reads.fastq.gz')
        with gzip.open(reads, 'wt') as output:
            for i in range(20):
                output.write('@duplicate\n{}\n{}\n+\n{}\n{}\n'.format(
                    sequence[:250], sequence[250:], 'I' * 250, 'I' * 250))
        more_reads = os.path.join(directory, 'more.fa.bz2')
        with bz2.open(more_reads, 'wt') as output:
            output.write('>duplicate\n{}\n{}\n'.format(sequence[:250], sequence[250:]))
        out = os.path.join(directory, 'assembly')
        command = [megahit, '-r', reads + ',' + more_reads, '-o', out, '-t', '1',
                   '-m', '800000000', '--mem-flag', '2',
                   '--k-list', '21', '--keep-tmp-files', '--overlap-correct']
        failed_env = env.copy()
        failed_env['MOCK_TRUNCATE'] = '1'
        run(command, failed_env, success=False, message='refusing partial correction')
        assert not os.path.exists(os.path.join(out, 'corrected.reads.fa'))
        run(command + ['--force'], env, success=False,
            message='.corrected.reads.partial.fa already exists')
        run([megahit, '--continue', '-o', out], env)
        corrected = os.path.join(out, 'corrected.reads.fa')
        with open(corrected) as source:
            lines = source.read().splitlines()
        assert len(lines) == 42
        assert lines[0] == '>mh0r' and lines[1] == sequence
        assert lines[-2] == '>mh20r' and lines[-1] == sequence
        with open(os.path.join(out, 'final.contigs.fa')) as source:
            assert source.read().startswith('>')
        run(command + ['--force'], env, success=False,
            message='corrected.reads.fa already exists')

        run([megahit, '-r', reads, '-o', os.path.join(directory, 'paired'),
             '--overlap-correct', '-1', reads, '-2', reads],
            env, success=False, message='requires single-end')
        missing_env = env.copy()
        missing_env['PATH'] = directory
        os.symlink(sys.executable, os.path.join(directory, 'python3'))
        run([megahit, '-r', reads, '-o', os.path.join(directory, 'missing'),
             '--overlap-correct'], missing_env, success=False,
            message='requires minimap2 on PATH')
        os.symlink('/bin/cat', os.path.join(directory, 'cat'))
        checkpoint = os.path.join(out, 'checkpoints.txt')
        with open(checkpoint, 'w') as output:
            output.write('0\tdone\n')
        run([megahit, '--continue', '-o', out], missing_env)
        malformed = os.path.join(directory, 'malformed.fq')
        with open(malformed, 'w') as output:
            output.write('@truncated\nACGT\n+\nII\n')
        run([megahit, '-r', malformed, '-o', os.path.join(directory, 'malformed'),
             '--overlap-correct'], env, success=False, message='Truncated FASTQ quality')
        no_overlaps_env = env.copy()
        no_overlaps_env['MOCK_NO_OVERLAPS'] = '1'
        run([megahit, '-r', reads, '-o', os.path.join(directory, 'no-overlaps'),
             '--overlap-correct'], no_overlaps_env, success=False,
            message='No read overlaps found')


if __name__ == '__main__':
    main()
