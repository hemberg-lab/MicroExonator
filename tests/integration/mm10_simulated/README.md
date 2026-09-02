# Lightweight mm10 smoke fixtures

These opt-in integration fixtures exercise full MicroExonator discovery and
quantification against the mouse mm10 reference. They are deliberately absent
from the default fast test suite's pipeline-execution path because a full mm10
FASTA, the complete GENCODE VM25 BED12 annotation, and the pipeline Conda
environments are required.

The two layouts share 100 truth events from the manuscript's Supplementary
Table 3: 65 VM25-annotated and 35 clean, strictly intronic non-annotated
microexons. This composition avoids treating annotated exon-boundary
extensions as novel cassette events. Each layout
has three groups (`A`, `B`, and `C`) and three biological replicates per group.
The single-end fixture is unstranded 100 nt; the paired fixture is FR 2x100 nt.
Both contain approximately 10% unrelated annotated exon-junction background.

## Contents

- `truth/events.tsv`: event coordinates, annotation status, flanking exons,
  profiles, target PSI values, and template hashes.
- `truth/*_reads.tsv.gz`: one row per single read or pair with isoform,
  anchors, insert size, orientation, substitution positions, and evidence
  validity.
- `single_end/fastq/`: nine single-end FASTQs.
- `paired_end/fastq/`: eighteen synchronized mate FASTQs.
- Per-layout `local_samples.tsv`, `bulk_samples.tsv`, and configuration
  templates; paired-end also has a headerless `paired_samples.tsv`.
- `PROVENANCE.json` and `SHA256SUMS`: source hashes, fixed seeds, and artifact
  integrity.

Coordinates are mm10, zero-based, half-open, matching BED and MicroExonator's
`chrom_strand_start_end` event representation. The complete manuscript table,
reference FASTA, VM25 annotation, and temporary generator are not committed.

## Fast integrity checks

The normal Python test suite validates gzip/FASTQ structure, 100 nt read
lengths, unique and synchronized identifiers, Phred ranges, event composition,
read-truth counts, anchor diversity, the >=12 valid included and skipped read
invariant for every event/sample, insert sizes, error/quality association, and
checksums:

```bash
python3 -m unittest tests.test_mm10_simulated_fixture -v
```

## Full reference-dependent smoke tests

Use Snakemake 7.32.4. Run single-end first and proceed to paired-end only after
it passes. Each command needs a new or empty work directory. Failed run
directories and `logs/snakemake.log` are retained for diagnosis.

```bash
python3 tests/integration/mm10_simulated/run_smoke_test.py \
  --layout single_end \
  --genome-fasta /path/to/mm10.fa \
  --annotation-bed /path/to/GENCODE_VM25_comprehensive.bed.gz \
  --snakemake /path/to/snakemake \
  --workdir /tmp/microexonator-mm10-single

python3 tests/integration/mm10_simulated/run_smoke_test.py \
  --layout paired_end \
  --genome-fasta /path/to/mm10.fa \
  --annotation-bed /path/to/GENCODE_VM25_comprehensive.bed.gz \
  --snakemake /path/to/snakemake \
  --workdir /tmp/microexonator-mm10-paired
```

The runner decompresses VM25 inside the work directory, builds a temporary
FASTA index there without modifying or copying the 2.6 GB reference, writes an
isolated concrete configuration, runs the pinned workflow, and invokes
`validate_smoke.py`. A passing run must report 99 of the 100 truth coordinates
in `Report/out.robustly_detected.txt` and produce numeric PSI, `CI_Lo`, and
`CI_Hi` for each reported truth event in all nine biological samples. The sole
permitted omission is `chr17_-_30598522_30598525` (ME053), a repeated 3-nt CAG
within one intron that the current reporting stage collapses into another CAG
coordinate. The validator rejects any other missing coordinate, including a
second missing truth event. Additional discoveries and PSI deviations are
informational in this first smoke-test version.
