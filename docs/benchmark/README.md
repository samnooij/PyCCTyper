# PyCCTyper benchmark

## Test details

To test computational performance of the original CCTyper against the
adapted PyCCTyper, we assembled two tests datasets:

- 100 random bacterial genomes with CRISPR-Cas arrays (`100_crisprs`)
  - this set contains 101 CRISPR-Cas loci
- 100 random genomes of _Campylobacter jejuni_ (n=75) and _C. coli_ (n=25; `100_campys`)
  - this set contains 56 CRISPR-Cas loci

Each dataset was analysed three times (3x) with each tool version,
using 1, 4 or 8 CPU threads. For each test, we recorded the total
runtime ('wallclock time'), mean CPU usage, maximum RAM use and
speed of ORF prediction (Prodigal or Pyrodigal, in minutes:seconds),
and _Cas_ gene identification (HMMer or PyHMMER; in minutes:seconds,
and iterations per second).

These metrics were recorded using the Linux shell command `/usr/bin/time -v`,
together with CCTyper's built-in logging mechanism which records the
time of separate processing steps in seconds.

Last but not least, output files from CCTyper, notably `CRISPR_Cas.tab`,
were compared between versions to verify identical outputs.

### System details

Benchmark tests were run on a HP EliteBook 840 G8 Notebook PC, with an
11th Generation Intel Core i7-1165G7 8-core processor and 16GB RAM.
As operating system, it uses Fedora Linux 44 (Workstation Edition)
with the GNOME desktop environment (version 50) and Linux kernel
7.1.6-201.fc44.x86_64.

### Software details

For the benchmark, we used CCTyper version 1.8.0, with dependencies
Prodigal v2.6.2 and HMMer v3.4. In the process of switching to 'PyCCTyper',
Prodigal was switched for Pyrodigal v3.7.1, and HMMer for PyHMMER v0.12.1.
For these two steps, we ran all benchmark tests separately.

## Results

### Identical outputs

CRISPR-Cas results were identical for the `100_crisprs` dataset.
For the `100_campy` dataset, however, there was a discrepancy:
56 CRISPR-Cas loci were consistently reported by all methods, but
a CRISPR array on contig `SAMN14992989.contig00028` was reported
3 out of 9 times: one of which by the regular CCTyper using 8 CPU threads,
and 2 by the version running Pyrodigal (1 and 4 threads).
This contig contains a CRISPR array with 6 repeats of length 36bp and 5 spacers
of length 30bp. The _Cas_ operon consists of only a Cas9_0_II gene.
Both the CRISPR array and _Cas_ operon (gene) are consistently identified
and reported in the respective `crisprs_all.tab` and `cas_operons_putative.tab`
files, however, 6/9 times the CRISPR subtype is predicted to be V-F1 with
a probability of 0.555, whereas in the 'correct three' it is II-C with 0.979
probability. This is likely where the problem arises: the prediction module
has some randomness built-in and a CRISPR-Cas locus is only reported
when both the CRISPR array (repeats) and _Cas_ operon are classified
as the same subtype.

So to conclude, in principle these different CCTyper variants produce
identical outputs. However, some randomness in the subtyping module
may cause inconsistent results, which may be associated with CRISPR
repeat sequences of subtype II-C that are a bit like subtype V-F1.

### Speed

### Memory use

## Details of ambiguous CRISPR subtyping result

The ambiguously identified/classified CRISPR-Cas locus had as repeat sequence:
`ATTTTACCATGTAAACAATTAATAATAGGCTAAAAC`. The online
[CRISPRCasdb](https://crisprcas.i2bc.paris-saclay.fr)
allows searching for repeat or spacer sequences using the BLAST algorithm,
which we used to investigate this sequence.
The sequence matches known repeats from _Campylobacter coli_ (3x) and
_Campylobacter jejuni_ (1x), that habe been classified as CRISPR-Cas
type II-C.
However, compared to the known repeat sequence, there is a mismatch at
position 11: a A->G mutation. The reference has an A, while the
query sequence has a G.
Searching the reference repeat sequence in CRISPRCasdb returns hits
with variants with 3 and 5 mismatches from _C. coli_, suggesting
that alternative variants may indeed occur and that repeat
sequences like this probably derive from type II-C CRISPR arrays
in _C. coli._

Looking further for this sequence in the AllTheBacteria dataset of
all _C. jejuni_ and _C. coli_ genomes reveals that this sequence is
found more often (n=154, _C. coli_ n=102, _C. jejuni_ n=51) and is always
reported as type V-F1 with confidence 0.555-0.737. The arrays have 4-11
repeats: a random mutation in 11 consecutive repeats seems unlikely, so this
is more likely a real repeat variant that has not been described before.
Also, 133 of these contigs have a _Cas9_ gene reported as type II-A/II-C
'orphan operon', suggesting that these form a functional but lesser known
CRISPR-Cas locus.
(The variant with an A is found 1156 times: ~7 times more often.
_C. coli_ n = 452, _C. jejuni_ n = 550.)

### The workflow in practice

To run the benchmarks and keep records of the commands and output,
I used the commands as listed below. In the git repository I
made different branches for the new features to implement (Pyrodigal and
PyHMMER). The workflow was as follows:

1. Check out git branch
2. Install current 'version' of CCTyper with `setup.py` file
3. Run test
4. Store and check results
5. Repeat!

Example code:

```bash
git checkout main
pip install .
/usr/bin/time -v cctyper test_data/100_campys.fasta results/100_campy-main-4\
 -t 4 --simplelog --db data | tee log/100_campy-main-4.log

git checkout pyrodigal
pip install .
/usr/bin/time -v cctyper test_data/100_campys.fasta results/100_campy-pyrodigal-4\
 -t 4 --simplelog --db data | tee log/100_campy-pyrodigal-4.log

pip checkout pyhmmer
pip install .
/usr/bin/time -v cctyper test_data/100_campys.fasta results/100_campy-pyhmmer-4\
 -t 4 --simplelog --db data | tee log/100_campy-pyhmmer-4.log
```

The results are saved under `docs/benchmark/`.
