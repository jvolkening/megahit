# Noisy viral assembly benchmark

This is a reproducible exploratory benchmark, not a model of all nanopore
errors or a validation of nanopore assembly. It simulates independent
substitution, insertion and deletion events at equal rates, with uniformly
distributed template lengths of 400-800 bases. Reads are sampled from linear
references, so coverage falls at their ends. The `--coverage` value is nominal:
the number of reads is `reference_length * coverage / 600`, rounded down.
There is a 15x random-sequence small-virus control and a 25x random-sequence
background control; the latter is not a realistic bacterial genome.

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
The benchmark creates hundreds of megabytes of reads and can take substantial
time and memory. Without `--real-references`, the generator uses random
sequences with approximately corresponding genome and segment lengths.

The evaluator uses sensitive minimap2 `map-ont` alignment and reports primary
alignments of at least 100 bases and 80% identity, reference coverage,
unaligned contigs, cross-reference alignments, and redundant aligned reference
bases. Redundancy is an overlap proxy, not a spurious-contig count: legitimate
repeats and overlapping contigs also contribute. Identity excludes unaligned
bases. Full misassembly, variant, false-positive and strain-specific analysis
requires additional tools and datasets.
