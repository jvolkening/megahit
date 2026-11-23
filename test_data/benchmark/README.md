# Noisy viral assembly benchmark

This is a reproducible exploratory benchmark, not a model of all nanopore
errors or a validation of nanopore assembly. It simulates independent
substitution, insertion and deletion events at equal rates, with uniformly
distributed template lengths of 400-800 bases. Reads are sampled from linear
references, so coverage falls at their ends. The `--coverage` value is nominal:
the number of reads is `reference_length * coverage / 600`, rounded down.
There is a 15x random-sequence small-virus control and a 25x random-sequence
background control by default; their depths can be changed with
`--virus-control-depth` and `--background-depth`. The background is not a
realistic bacterial genome.

The real-reference mode uses influenza A/California/07/2009 segments
NC_026431.1-NC_026438.1, SARS-CoV-2 NC_045512.2, and HSV-1 NC_001806.2.
Identical duplicated influenza records in the source are ignored. Four `R`
ambiguities in the influenza references are resolved to `A`. This model has no
systematic, homopolymer-associated or correlated nanopore errors. Optionally,
`--minor-strain-depth 60` adds a 2%-substitution-divergent copy of the first
influenza segment at nominal 60x coverage; this simple strain model does not
capture real diversity. Do not infer false-positive rates or
strain preservation on real datasets from these results.

To reproduce the reference inputs in an empty working directory, using `gh`:

```sh
gh api -H 'Accept: application/vnd.github.raw+json' 'repos/ncbi/ngs-tools/contents/tools/tax/examples/sequence_tree/Viruses/Influenzavirus%20A/California/07/2009%28H1N1%29%29/641809.fasta?ref=f8ed491f6635034fe8f48310ad6cc5da29c08ea5' > flu_source.fa
gh api -H 'Accept: application/vnd.github.raw+json' 'repos/onecodex/sars-cov-2/contents/reference/NC_045512.2.reference.fasta?ref=2fd8baffdd8ef71d2383d621cc1f78e43807f94e' > corona_source.fa
gh api -H 'Accept: application/vnd.github.raw+json' 'repos/yyoshiaki/VIRTUS/contents/data/NC_001806.2.fasta?ref=49faf55f97c8f3084b426d2db6640519d6f2ce71' > herpes_source.fa
```

With `MEGAHIT_SRC` set to the source directory, a built `megahit` on PATH and
minimap2 installed:

```sh
python3 "$MEGAHIT_SRC/test_data/benchmark/simulate.py" trial --real-references --coverage 1100 --error 0.10
megahit -r trial.reads.fa -o trial-out -t 2 -m 1800000000 --mem-flag 2
python3 "$MEGAHIT_SRC/test_data/benchmark/evaluate.py" minimap2 trial.ref.fa trial-out/final.contigs.fa
```

Repeat with `--error 0.05` and a different prefix for the lower-error control.
For a smaller related-strain experiment, add `--flu-only --minor-strain-depth 60`.
Use `--seed` with a different value and a new prefix for independent replicates.
To compare the opt-in local filter on the same reads, rerun MEGAHIT in a
different output directory with `--min-local-support 2`. The filter can remove
genuine minor strains, even when excluded candidates are retained separately;
it is not a strain-aware consensus algorithm.
In this mode, single-end reads are mapped once across threads during local
assembly; legacy mode can process each read once per thread. Support counts
are k-mer observations, not independent-read counts, and include the injected
parent contig.
In one 5%-error influenza simulation, filtering reduced 947 contigs to 38,
but a 60x minor strain's co-occurring SNPs were no longer detected in either
the final assembly or the excluded candidates. At 10% error the filter did not
reduce output in the corresponding experiment (139 to 311 contigs).
In a second 5%-error seed, filtering reduced 1,863 contigs to 43 and reduced
15x control coverage from 2,370 to 2,203 of 2,500 bases; neither assembly had
qualifying minor-strain evidence. These simulations use independent errors
and do not establish performance on real nanopore data.
Candidate identifiers may recur across iteration files; rename them before
combining those files for evaluation.

On the 5%-error mixed-genome dataset with two threads, the default assembly
produced 67,014 contigs and the experimental mode produced 59. However, HSV-1
reference coverage fell from 152,222 to 135,066 of 152,222 bases, and coverage
of the 15x control fell from 2,312 to 1,793 of 2,500 bases. The experimental
run also changes single-end read scheduling as described above; it used the
portable binary while the default used the BMI2 binary. These observations
do not establish strain-safe filtering or real-data performance.
At 10% error, the mixed-genome experimental run used one thread and produced
1,126 contigs, including 34 without a passing alignment; reference coverage
was 21,829/29,903 bases for SARS-CoV-2 and 82,801/152,222 for HSV-1. The
two-thread default run is not a controlled comparison because legacy local
read processing depends on thread count.

For a nominal 10x mixed community, including both controls, use:

```sh
python3 "$MEGAHIT_SRC/test_data/benchmark/simulate.py" trial10 --real-references --coverage 10 --virus-control-depth 10 --background-depth 10 --error 0.02
megahit -r trial10.reads.fa -o trial10-default -t 1 -m 800000000 --mem-flag 2
megahit -r trial10.reads.fa -o trial10-filtered -t 1 -m 800000000 --mem-flag 2 --min-local-support 2
python3 "$MEGAHIT_SRC/test_data/benchmark/evaluate.py" minimap2 trial10.ref.fa trial10-filtered/final.contigs.fa
```

Repeat with `--error 0.05` and a new prefix. In single-seed comparisons,
filtering reduced contigs from 298 to 49 at 2% error and from 651 to 356 at
5% error. It also reduced HSV-1 reference coverage from 139,993 to 135,664
of 152,222 bases at 2%, and from 144,014 to 122,633 at 5%. SARS-CoV-2
coverage fell from 29,459 to 29,301 of 29,903 bases at 2%, and from 29,394 to
25,830 at 5%. The simulator samples linear references, leaving their ends
underrepresented; these results do not establish real-data performance.

Filtered assembly by reference (summed lengths include overlaps and
unaligned ends, and are not unique reference coverage):

| Reference | Expected bp | 2% contigs | 2% bp | 5% contigs | 5% bp |
| --- | ---: | ---: | ---: | ---: | ---: |
| Influenza NC_026431.1 | 982 | 1 | 854 | 1 | 652 |
| Influenza NC_026432.1 | 863 | 1 | 773 | 1 | 679 |
| Influenza NC_026433.1 | 1,701 | 1 | 1,501 | 3 | 1,554 |
| Influenza NC_026434.1 | 1,410 | 1 | 1,070 | 4 | 1,346 |
| Influenza NC_026435.1 | 2,274 | 1 | 2,131 | 4 | 2,081 |
| Influenza NC_026436.1 | 1,497 | 1 | 1,283 | 2 | 963 |
| Influenza NC_026437.1 | 2,151 | 2 | 2,001 | 3 | 1,788 |
| Influenza NC_026438.1 | 2,280 | 1 | 2,103 | 2 | 1,709 |
| SARS-CoV-2 | 29,903 | 8 | 29,359 | 57 | 28,042 |
| HSV-1 | 152,222 | 26 | 135,600 | 228 | 130,324 |
| Synthetic virus | 2,500 | 1 | 2,195 | 4 | 1,828 |
| Synthetic bacterial background | 30,000 | 5 | 29,996 | 47 | 28,441 |
| No passing alignment | - | 0 | 0 | 0 | 0 |
| **Total** | - | **49** | **208,866** | **356** | **199,407** |

### Exploratory overlap correction at 10x

These are research comparisons, not recommended production settings. The
simulator's independent errors are easier than realistic nanopore errors, and
the results below expose potential strain loss and misassembly.

At 5% error, three seeds (932414, 932415, 932416) were assembled with one
thread and `--min-local-support 2`. Using the existing k-list beginning with
15 instead of 21 produced 255, 278 and 263 contigs, compared with 356, 298
and 347 for the original list. For the first seed, one 331 bp contig had two
passing primary alignments to distant HSV-1 positions on opposite strands
(a possible misjoin); cross-reference alignment metrics alone do not catch
such within-reference events. At 2% error, starting at k=15 produced 32 HSV-1
contigs versus 26 with k=21, so the change is not uniformly better.

As a separate experiment, one round of minimap2 v2.31 overlap mapping and
Racon v1.5.0 fragment consensus was applied before assembly with k=15 and
the local filter. The commands below reproduce the first seed from `trial10`
generated above with `--error 0.05`, using the same preprocessing parameters.
MEGAHIT's experimental
`--overlap-correct` option now performs these two preprocessing steps using
tools on `PATH`; it accepts single-end input only, normalizes identifiers,
discards FASTQ quality scores, and refuses incomplete corrected read sets.
It preserves the original files and writes `corrected.reads.fa` separately.
It does not set the k-list or local filter automatically. Renaming reads can
change overlap tie-breaking, so results need not exactly match the earlier
manual run; the integrated first-seed run produced 79 contigs and covered
135,199/152,222 HSV-1 bases.

```sh
megahit -r trial10.reads.fa -o trial10-corrected -t 1 -m 800000000 --mem-flag 2 --overlap-correct --min-local-support 2 --k-list 15,21,29,39,59,79,99,119,141
python3 "$MEGAHIT_SRC/test_data/benchmark/evaluate.py" minimap2 trial10.ref.fa trial10-corrected/final.contigs.fa
```

In the three 5% seeds, corrected-read assemblies produced 77, 73 and 76
contigs (N50 5,495, 6,049 and 5,875 bp), compared with 356, 298 and 347
contigs (N50 646, 688 and 663 bp) for uncorrected k=21 filtered assemblies.
For the first seed, corrected-read HSV-1 coverage was 135,304/152,222 bp
versus 122,633 bp in the uncorrected filtered assembly. The unfiltered,
uncorrected k=21 assembly covered 144,014 bp. Racon emitted 3,791 reads
for 3,791 input reads in each of these runs, but this is not guaranteed by
these settings for other inputs.

The method is **not strain-safe**: in a separate 10x major / 5x minor,
2%-divergent influenza simulation, eight original reads contained at least
four co-occurring minor SNPs (covering 19 of 20 sites); no corrected reads
met that heuristic. The uncorrected k=15 assembly without local filtering
had two qualifying contigs covering 16 sites, but neither a filtered
corrected assembly nor an assembly combining raw and corrected reads without
the filter had qualifying contigs. Those observations do not prove every
minor allele was erased, but they demonstrate loss of the tested haplotype
signal. As another exploratory tradeoff, selecting one raw or corrected
copy per read by a fixed SHA-256 hash retained no qualifying minor contigs
with 50% raw reads; 75% raw produced one qualifying contig with 15 sites
without filtering, but 107 contigs overall versus 83 for the all-raw
unfiltered run. Filtering lost the qualifying signal at both proportions.
These single-mixture results do not establish a safe correction policy.

Across the mixed-genome tests, 53-171 overlap pairs per replicate had
high-mapping-quality read alignments at disjoint positions on the same
reference; repeats or ambiguous mappings can produce these apparent
discordances. Reference alignment and overlap counts do not exclude
misassemblies. Keep the original reads, and do not treat consensus correction
or the local filter as preserving minority strains. All-vs-all mapping can
also become costly at high depth.

The benchmark creates hundreds of megabytes of reads and can take substantial
time and memory. Without `--real-references`, the generator uses random
sequences with approximately corresponding genome and segment lengths.

The evaluator uses sensitive minimap2 `map-ont` alignment and reports primary
alignments of at least 100 bases and 80% identity, reference coverage,
unaligned contigs, cross-reference alignments, and redundant aligned reference
bases. Per reference, it also reports the number and summed full lengths of
distinct assigned contigs, including their unaligned ends. A contig with
passing alignments to multiple references is counted under each.
Redundancy is an overlap proxy, not a spurious-contig count: legitimate
repeats and overlapping contigs also contribute. Identity excludes unaligned
bases. Full misassembly, variant, false-positive and strain-specific analysis
requires additional tools and datasets.

For the optional minor-strain simulation, `strain_evaluate.py` checks whether
at least four simulated minor SNPs occur together in the same primary
alignment. This is a heuristic specific to the substitution-only simulated
strain, not a general haplotype or false-positive evaluation:

```sh
python3 "$MEGAHIT_SRC/test_data/benchmark/strain_evaluate.py" minimap2 trial.ref.fa trial-out/final.contigs.fa NC_026431.1 NC_026431.1_minor
```
