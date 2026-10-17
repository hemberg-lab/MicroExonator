"""PSI-Sigma module: reference, per-run junction and intron-retention caches, comparisons.

PSI-Sigma v2.3 (Lin & Krainer, Bioinformatics 2019) was written for STAR BAMs.
The umbrella feeds it the shared HISAT2 alignment instead, with three
adaptations, each faithful to PSI-Sigma's own definitions:

  * uniqueness: PSI-Sigma treats MAPQ >= 255 (STAR) as unique; HISAT2's unique
    alignments carry NH:i:1 (MAPQ 60). Junction records are rewritten in the
    stream (MAPQ 255 when NH:i:1, else 0) before PSI-Sigma's own, unmodified
    awk counts them; intron coverage uses NH == 1 directly. No BAM is written.
  * mates: the awk counts a junction once per read name only when records of a
    name are adjacent, so the spliced records are grouped by name
    (`samtools collate`), as PSI-Sigma's driver does with `samtools sort -n`.
  * canonical junctions (on by default): PSI-Sigma's README aligns with STAR
    --outFilterIntronMotifs RemoveNoncanonical. Junctions whose motif is not
    GT-AG, GC-AG or AT-AC (either strand) are dropped from the junction cache.

Intron retention (--irmode 0, PSI-Sigma's default and the only mode allowed
here) depends only on the annotation: each retained-intron event's five
coverage positions come from the intron's own coordinates. So the per-run cache
keeps the coverage at every such position, counted with PSI-Sigma's read rules
(primary, proper pair for paired-end, unique, no N/S/H in the CIGAR; M bases
count, N and D advance). A comparison sums a replicate's runs, takes PSI-Sigma's
median per intron and never needs a BAM.

reference  annotation for PSI-Sigma (primary chromosomes; gene and transcript
           lines added where missing; sorted with transcript before exon), its
           transcript mapping, the intron-retention universe from PSI-Sigma's own
           database builder run without junctions, and a genome index
prep       per run: run.SJ.tab.gz (PSI-Sigma's 7 columns), ir_positions.tsv.gz,
           qc.json, cache.json
compare    per comparison: replicate sums, PSI-Sigma's database, PSI and filter
           steps (native, unmodified), normalized results (effect = PSI(A) -
           PSI(B) as a fraction), native.tar.gz, status.json
"""

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

if __package__:
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from umbrella_modules import NORMALIZED, write_status

PRIMARY = re.compile(r"^[A-Za-z0-9]+$")   # PSI-Sigma splits IDs on "_"; "." breaks its IR step
CANONICAL = {("GT", "AG"), ("GC", "AG"), ("AT", "AC"),     # + strand
             ("CT", "AC"), ("CT", "GC"), ("GT", "AT")}     # the same motifs read on - strand
NH_TO_MAPQ = r"""BEGIN{OFS="\t"} {u=0; for(i=12;i<=NF;i++) if($i=="NH:i:1"){u=1; break}; $5=(u?255:0); print}"""
NATIVE = ("PSIsigma-db-v.1.0.pl", "PSIsigma-PSI-v.1.1.pl", "PSIsigma-filter-v.1.0.pl",
          "sjFromSAMcollapseUandM_inclOverlaps.awk")
EFFECT = "PSI-Sigma dPSI (%) / 100: mean PSI(A) - mean PSI(B)"
STATISTIC = "Welch t-test (PSI-Sigma --variance 0 runs PDL::Stats t_test_nev); BH FDR over all events"
FILTER_HEADER = ("Event Region", "Gene Symbol", "Target Exon", "Event Type", "N", "T", "Exon Type",
                 "Reference Transcript", "dPSI (%)", "T-test p-value", "FDR (BH)", "N Values",
                 "T Values", "Database ID")
TABLES = ("results.tsv.gz", "native.tar.gz")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command, log, cwd=None, env=None, shell=False):
    log.write("$ " + (command if shell else " ".join(map(str, command))) + "\n")
    log.flush()
    subprocess.run(command if shell else [str(c) for c in command], check=True, cwd=cwd, env=env,
                   stdout=log, stderr=subprocess.STDOUT, shell=shell,
                   executable="/bin/bash" if shell else None)


# ------------------------------------------------------------ reference

def _attribute(attributes, key):
    match = re.search(r'(?:^|;\s*){} "([^"]*)"'.format(key), attributes)
    return match.group(1) if match else None


def prepare_gtf(source, target):
    """PSI-Sigma's annotation: primary chromosomes; every transcript with a
    transcript line and every gene with a gene line (spans from its exons when
    the source has none); sorted by chromosome and start with gene, then
    transcript, then exon at equal starts (PSI-Sigma's database builder walks
    transcript lines to group loci; an exon sorted before its transcript would
    be given to the previous locus)."""
    genes, transcripts, exons, kept, dropped = {}, {}, [], set(), set()
    with _open(source) as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("gene", "transcript", "exon"):
                continue
            if not PRIMARY.match(fields[0]):
                dropped.add(fields[0])
                continue
            kept.add(fields[0])
            start, end = int(fields[3]), int(fields[4])
            gene_id = _attribute(fields[8], "gene_id")
            if fields[2] == "gene":
                genes[gene_id] = fields
                continue
            transcript_id = _attribute(fields[8], "transcript_id")
            if not gene_id or not transcript_id:
                raise SystemExit("PSI-Sigma annotation line without gene_id or transcript_id: " + line[:200])
            if fields[2] == "transcript":
                transcripts[transcript_id] = fields
                continue
            exons.append(fields)
            spans = transcripts.setdefault("~" + transcript_id, [fields[0], start, end, fields[6], gene_id,
                                                                 fields[8]])
            spans[1], spans[2] = min(spans[1], start), max(spans[2], end)
    lines = []
    synthesized = {"transcript": 0, "gene": 0}
    gene_spans = {}
    for key, value in transcripts.items():
        if key.startswith("~"):
            transcript_id = key[1:]
            chrom, start, end, strand, gene_id, attributes = value
            if transcript_id not in transcripts:
                gene_name = _attribute(attributes, "gene_name")
                fields = [chrom, "umbrella", "transcript", str(start), str(end), ".", strand, ".",
                          'gene_id "{}"; transcript_id "{}";{}'.format(
                              gene_id, transcript_id, ' gene_name "{}";'.format(gene_name) if gene_name else "")]
                for attribute in ("transcript_type", "transcript_biotype"):
                    value_ = _attribute(attributes, attribute)
                    if value_:
                        fields[8] += ' {} "{}";'.format(attribute, value_)
                lines.append(fields)
                synthesized["transcript"] += 1
            span = gene_spans.setdefault(gene_id, [chrom, start, end, strand, attributes])
            span[1], span[2] = min(span[1], start), max(span[2], end)
        else:
            lines.append(value)
    for gene_id, fields in genes.items():
        lines.append(fields)
    for gene_id, (chrom, start, end, strand, attributes) in gene_spans.items():
        if gene_id in genes:
            continue
        gene_name = _attribute(attributes, "gene_name")
        lines.append([chrom, "umbrella", "gene", str(start), str(end), ".", strand, ".",
                      'gene_id "{}";{}'.format(gene_id, ' gene_name "{}";'.format(gene_name) if gene_name else "")])
        synthesized["gene"] += 1
    lines.extend(exons)
    order = {"gene": 0, "transcript": 1, "exon": 2}
    lines.sort(key=lambda f: (f[0], int(f[3]), order[f[2]], int(f[4]), f[8]))
    with open(target, "w") as stream:
        for fields in lines:
            stream.write("\t".join(fields) + "\n")
    return {"chromosomes": sorted(kept), "dropped_contigs": len(dropped), "synthesized": synthesized,
            "exons": len(exons)}


def write_mapping(gtf, target):
    """PSI-Sigma's transcript -> gene name, strand file, as its driver writes it."""
    rows = 0
    with open(gtf) as stream, open(target, "w") as out:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "transcript":
                continue
            transcript = _attribute(fields[8], "transcript_id").replace("_", ".")
            name = (_attribute(fields[8], "gene_name") or _attribute(fields[8], "gene_id")).replace("_", ".")
            out.write("{}\t{}\t{}\n".format(transcript, name, fields[6]))
            rows += 1
    if not rows:
        raise SystemExit("PSI-Sigma annotation has no transcript lines")


def gtf_chromosomes(gtf):
    chromosomes = []
    with open(gtf) as stream:
        for line in stream:
            chrom = line.split("\t", 1)[0]
            if chrom not in chromosomes:
                chromosomes.append(chrom)
    return [c for c in chromosomes if not re.match(r"^(chr)?(GL|KI)", c)]


def build_database(psisigma, gtf, chromosomes, groupa, groupb, name, cwd, threads, log, env=None):
    """PSI-Sigma's database builder, per chromosome (novel junctions >= 10
    unique reads, short reads, irmode 0), concatenated as its driver does."""
    def one(chrom):
        with open(Path(cwd) / ("db.{}.log".format(chrom)), "w") as chrom_log:
            run(["perl", Path(psisigma) / "PSIsigma-db-v.1.0.pl", Path(gtf).resolve(), chrom, 10, 1, 0,
                 groupa, groupb], chrom_log, cwd=cwd, env=env)
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        list(pool.map(one, chromosomes))
    database = Path(cwd) / (name + ".db")
    with open(database, "w") as out:
        for chrom in chromosomes:
            part = Path(cwd) / "{}.{}.db.tmp".format(name, chrom)
            if not part.exists():
                raise SystemExit("PSI-Sigma database builder wrote no {}".format(part.name))
            with open(part) as stream:
                shutil.copyfileobj(stream, out)
    with open(database) as stream:
        log.write("database: {} lines\n".format(sum(1 for _ in stream)))
    return database


def ir_positions(database):
    """{(chrom, i2s, i2e): sorted positions} from the R (retained intron) lines
    of a PSI-Sigma database, exactly as PSIsigma-ir-v.1.2.pl derives them."""
    universe = {}
    with open(database) as stream:
        for line in stream:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, i2s, i2e, start, end, name = fields[0], fields[3], fields[4], int(fields[8]), int(fields[9]), fields[10]
            if "_R_" not in name:
                continue
            interval = int((end - start) / 5)
            positions = {start + 1, end - 1, start + interval - 1, start + 2 * interval - 1,
                         start + 3 * interval - 1}
            universe.setdefault((chrom, i2s, i2e), set()).update(positions)
    return {key: sorted(value) for key, value in universe.items()}


def write_universe(universe, path):
    with gzip.open(path, "wt") as stream:
        for (chrom, i2s, i2e), positions in sorted(universe.items(), key=lambda kv: (kv[0][0], int(kv[0][1]),
                                                                                       int(kv[0][2]))):
            stream.write("{}\t{}\t{}\t{}\n".format(chrom, i2s, i2e, ",".join(map(str, positions))))


def read_universe(path):
    universe = {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            chrom, i2s, i2e, positions = line.rstrip("\n").split("\t")
            universe[(chrom, i2s, i2e)] = [int(p) for p in positions.split(",")]
    return universe


def reference(gtf, genome, psisigma, out_dir, threads=1, record=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    prepared = out / "annotation.psisigma.gtf"
    summary = prepare_gtf(gtf, prepared)
    write_mapping(prepared, out / "annotation.psisigma.gtf.mapping.txt")
    with open(out / "reference.log", "w") as log, tempfile.TemporaryDirectory(prefix="universe-", dir=out) as directory:
        work = Path(directory)
        # two samples without junctions: the database then holds annotation events only
        for label in ("u1", "u2"):
            (work / (label + ".SJ.out.tab")).write_text("")
        (work / "universe.groupa.txt").write_text("u1\n")
        (work / "universe.groupb.txt").write_text("u2\n")
        database = build_database(psisigma, prepared, gtf_chromosomes(prepared), "universe.groupa.txt",
                                  "universe.groupb.txt", "universe", work, threads, log)
        universe = ir_positions(database)
    if not universe:
        raise SystemExit("PSI-Sigma found no retained-intron events in the annotation")
    write_universe(universe, out / "ir_universe.tsv.gz")
    fai = out / "genome.fa.fai"
    subprocess.run(["samtools", "faidx", str(genome), "--fai-idx", str(fai)], check=True)
    Path(record or out / "reference.json").write_text(json.dumps(dict(
        summary, introns=len(universe), positions=sum(len(p) for p in universe.values()),
        annotation_sha256=_sha256(prepared), universe_sha256=_sha256(out / "ir_universe.tsv.gz")),
        sort_keys=True, indent=2) + "\n")


# ------------------------------------------------------------ per run: junctions

def canonical(sequence_of, chrom, start, end):
    """Motif of an intron (1-based closed): its first two and last two bases."""
    donor = sequence_of(chrom, start - 1, start + 1).upper()
    acceptor = sequence_of(chrom, end - 2, end).upper()
    return (donor, acceptor) in CANONICAL


def junctions(bam, awk, out, tmp, threads, log):
    """PSI-Sigma's 7-column junction counts from the spliced records, grouped by
    name, with HISAT2 uniqueness (NH:i:1) given as MAPQ 255."""
    prefix = Path(tmp) / "collate"
    command = ("set -o pipefail; samtools view -h -F 4 {bam} "
               "| awk -F'\\t' '/^@/ || $6 ~ /N/' "
               "| samtools collate -O -u -@ {threads} -T {prefix} - "
               "| samtools view - | awk -F'\\t' '{nh}' | awk -f {awk} > {out}").format(
        bam=_quote(bam), threads=max(1, threads - 1), prefix=_quote(prefix), nh=NH_TO_MAPQ,
        awk=_quote(awk), out=_quote(out))
    run(command, log, shell=True)


def _quote(path):
    return "'" + str(path).replace("'", "'\\''") + "'"


def filter_canonical(raw, target, fasta=None, fai=None, enabled=True):
    """Sorted, gzipped junction cache; non-canonical junctions dropped when enabled."""
    rows = []
    with open(raw) as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 7:
                rows.append(fields)
    removed = {"junctions": 0, "unique_fragments": 0}
    if enabled and rows:
        import pysam
        genome = pysam.FastaFile(str(fasta), filepath_index=str(fai)) if fai else pysam.FastaFile(str(fasta))
        kept = []
        for fields in rows:
            if canonical(lambda c, s, e: genome.fetch(c, s, e), fields[0], int(fields[1]), int(fields[2])):
                kept.append(fields)
            else:
                removed["junctions"] += 1
                removed["unique_fragments"] += int(fields[3])
        rows = kept
    rows.sort(key=lambda f: (f[0], int(f[1]), int(f[2])))
    with gzip.open(target, "wt") as stream:
        for fields in rows:
            stream.write("\t".join(fields) + "\n")
    return dict(removed, kept=len(rows))


# ------------------------------------------------------------ per run: intron coverage

def qualifies(flag, nh, cigar, paired):
    """PSI-Sigma's intron-read rule (PSIsigma-ir-v.1.2.pl, short reads):
    `samtools view -F 256 [-f 2] -q 255`, then no N, S or H in the CIGAR; with
    HISAT2, -q 255 is NH:i:1. Unmapped records are never fetched by region."""
    if flag & 0x4 or flag & 0x100:
        return False
    if paired and not flag & 0x2:
        return False
    if nh != 1:
        return False
    return not any(op in (3, 4, 5) for op, _ in cigar)


def blocks(position, cigar):
    """1-based half-open reference intervals PSI-Sigma counts: M bases count;
    N and D advance without counting; anything else (I, P, and, as in
    PSI-Sigma, = and X) neither counts nor advances."""
    start = position
    for op, length in cigar:
        if op == 0:
            yield start, start + length
            start += length
        elif op in (2, 3):
            start += length


class PositionCounter:
    """Coverage at a fixed set of positions per chromosome, by difference
    arrays over the sorted positions (no per-base work)."""

    def __init__(self, universe):
        self.positions = {}
        for (chrom, _, _), positions in universe.items():
            self.positions.setdefault(chrom, set()).update(positions)
        self.positions = {chrom: sorted(p) for chrom, p in self.positions.items()}
        self.diff = {chrom: [0] * (len(p) + 1) for chrom, p in self.positions.items()}

    def add(self, chrom, intervals):
        positions = self.positions.get(chrom)
        if not positions:
            return
        diff = self.diff[chrom]
        for start, end in intervals:
            low = bisect.bisect_left(positions, start)
            high = bisect.bisect_left(positions, end)
            if low < high:
                diff[low] += 1
                diff[high] -= 1

    def counts(self):
        result = {}
        for chrom, positions in self.positions.items():
            running, diff = 0, self.diff[chrom]
            for index, position in enumerate(positions):
                running += diff[index]
                if running:
                    result[(chrom, position)] = running
        return result


def coverage(bam, universe, paired, qc):
    """Count PSI-Sigma's intron reads at the universe positions, and the QC
    counts that show whether HISAT2's flags behave as the adapter assumes."""
    import pysam
    counter = PositionCounter(universe)
    with pysam.AlignmentFile(str(bam)) as alignments:
        names = set(alignments.references)
        for chrom in sorted(counter.positions):
            if chrom not in names:
                continue
            for read in alignments.fetch(chrom):
                flag = read.flag
                qc["records"] += 1
                nh = read.get_tag("NH") if read.has_tag("NH") else None
                if nh is None:
                    qc["missing_nh"] += 1
                cigar = read.cigartuples or []
                primary = not flag & 0x904
                if primary:
                    qc["primary"] += 1
                    if flag & 0x2:
                        qc["proper_pair"] += 1
                    if any(op == 4 for op, _ in cigar):
                        qc["primary_soft_clipped"] += 1
                if read.mapping_quality == 60 and nh not in (None, 1):
                    qc["mapq60_not_unique"] += 1
                if nh == 1 and read.mapping_quality != 60:
                    qc["unique_not_mapq60"] += 1
                if any(op in (7, 8) for op, _ in cigar):
                    qc["eq_x_cigar"] += 1
                if qualifies(flag, nh, cigar, paired):
                    qc["intron_reads"] += 1
                    counter.add(chrom, blocks(read.reference_start + 1, cigar))
    return counter.counts()


def prep(bam, universe_path, psisigma, out_dir, run_id, layout, fasta, fai=None, threads=2,
         canonical_only=True, identity=None, tmp_root=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if layout not in ("PE", "SE"):
        raise SystemExit("layout must be PE or SE")
    universe = read_universe(universe_path)
    qc = dict.fromkeys(("records", "primary", "proper_pair", "primary_soft_clipped", "missing_nh",
                        "mapq60_not_unique", "unique_not_mapq60", "eq_x_cigar", "intron_reads"), 0)
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out)
    with tempfile.TemporaryDirectory(prefix="psisigma-prep-", dir=tmp_root) as directory, \
            open(out / "prep.log", "w") as log:
        work = Path(directory)
        with ThreadPoolExecutor(max_workers=2) as pool:
            raw = work / "run.SJ.out.tab"
            sj = pool.submit(junctions, bam, Path(psisigma) / "sjFromSAMcollapseUandM_inclOverlaps.awk",
                             raw, work, threads, log)
            counts = coverage(bam, universe, layout == "PE", qc)
            sj.result()
        filtered = filter_canonical(raw, work / "run.SJ.tab.gz", fasta, fai, canonical_only)
        if qc["missing_nh"]:
            raise SystemExit("{} alignments of {} have no NH tag; HISAT2 uniqueness cannot be "
                             "given to PSI-Sigma".format(qc["missing_nh"], run_id))
        if layout == "PE" and qc["primary"] and not qc["proper_pair"]:
            raise SystemExit("{} is paired-end in the manifest but has no proper pairs".format(run_id))
        if layout == "SE" and qc["proper_pair"]:
            raise SystemExit("{} is single-end in the manifest but has proper pairs".format(run_id))
        with gzip.open(work / "ir_positions.tsv.gz", "wt") as stream:
            for (chrom, position), count in sorted(counts.items()):
                stream.write("{}\t{}\t{}\n".format(chrom, position, count))
        shutil.move(str(work / "run.SJ.tab.gz"), str(out / "run.SJ.tab.gz"))
        shutil.move(str(work / "ir_positions.tsv.gz"), str(out / "ir_positions.tsv.gz"))
    qc["soft_clip_fraction"] = round(qc["primary_soft_clipped"] / qc["primary"], 6) if qc["primary"] else 0.0
    qc["junctions"] = filtered
    (out / "qc.json").write_text(json.dumps(qc, sort_keys=True, indent=2) + "\n")
    (out / "cache.json").write_text(json.dumps(dict(
        identity or {}, tool="psisigma", run_id=run_id, layout=layout, canonical_only=bool(canonical_only),
        universe_sha256=_sha256(universe_path),
        sha256={name: _sha256(out / name) for name in ("run.SJ.tab.gz", "ir_positions.tsv.gz")}),
        sort_keys=True, indent=2) + "\n")


# ------------------------------------------------------------ per comparison

def median_value(values):
    """PSIsigma-ir-v.1.2.pl median, then sprintf("%d"): truncation toward zero."""
    values = sorted(values)
    middle = len(values) // 2
    value = values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2
    return int(value)


def replicate_tables(runs, caches, universe):
    """Junction and intron tables of one biological replicate: summed over its
    runs, as a BAM merged from them would give."""
    junction_sums, position_sums = {}, {}
    for run in runs:
        with gzip.open(Path(caches[run]) / "run.SJ.tab.gz", "rt") as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                key = tuple(fields[:3])
                current = junction_sums.setdefault(key, [0, 0, 0, 0])
                for index in range(4):
                    current[index] += int(fields[3 + index])
        with gzip.open(Path(caches[run]) / "ir_positions.tsv.gz", "rt") as stream:
            for line in stream:
                chrom, position, count = line.rstrip("\n").split("\t")
                key = (chrom, int(position))
                position_sums[key] = position_sums.get(key, 0) + int(count)
    junction_lines = ["\t".join(list(key) + [str(v) for v in values]) + "\n"
                      for key, values in sorted(junction_sums.items(), key=lambda kv: (kv[0][0], int(kv[0][1]),
                                                                                         int(kv[0][2])))]
    intron_lines = ["{}\t{}\t{}\t{}\n".format(chrom, i2s, i2e, median_value(
        [position_sums.get((chrom, p), 0) for p in positions]))
        for (chrom, i2s, i2e), positions in sorted(universe.items())]
    return junction_lines, intron_lines


def check_cache(record, cache, run, ids, settings, universe_sha256):
    expected = {"reference_id": ids.get("reference_id"), "module_reference_id": ids.get("module_reference_id"),
                "cache_id": Path(cache).name if ids else None, "universe_sha256": universe_sha256,
                "canonical_only": settings.get("canonical_only")}
    for key, value in expected.items():
        if value is not None and record.get(key) != value:
            raise ValueError("PSI-Sigma cache {} of {} has {} {}, expected {}".format(
                cache, run, key, record.get(key), value))
    for name, digest in record.get("sha256", {}).items():
        if _sha256(Path(cache) / name) != digest:
            raise ValueError("PSI-Sigma cache {} of {}: checksum mismatch for {}".format(cache, run, name))


def compare(preflight, caches, reference_dir, psisigma, out_dir, threads=1, ids=None, settings=None,
            perl5lib=None):
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    settings = dict({"nread": 10, "skipratio": 0.05, "irrange": 5, "variance": 0, "irclean": 0,
                     "canonical_only": True}, **(settings or {}))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ids = ids or {}
    base = dict(ids, project_id=preflight.get("project_id", ""), comparison_id=preflight.get("comparison_id", ""))
    reasons = []
    if not preflight.get("inference_supported", False):
        reasons += preflight.get("reasons", []) or ["unsupported comparison design"]
    if any(len(preflight["replicates"][side]) < 2 for side in ("a", "b")):
        reasons.append("PSI-Sigma's t-test needs at least two biological replicates per side")
    if reasons:
        _empty_outputs(out)
        write_status(out / "status.json", "psisigma", "unsupported", reasons, **base)
        return
    reference_dir = Path(reference_dir)
    universe_path = reference_dir / "ir_universe.tsv.gz"
    universe = read_universe(universe_path)
    universe_sha256 = _sha256(universe_path)
    members = {}
    for run_, replicate in preflight["collapse"].items():
        members.setdefault(replicate, []).append(run_)
    for run_ in preflight["collapse"]:
        record = json.loads((Path(caches[run_]) / "cache.json").read_text())
        if record.get("run_id") != run_:
            raise ValueError("PSI-Sigma cache {} belongs to {}, not {}".format(caches[run_], record.get("run_id"), run_))
        check_cache(record, caches[run_], run_, ids, settings, universe_sha256)
    # PSI-Sigma splits sample names on "_" and "."; plain labels avoid that.
    # Its groupa is N and groupb T, and dPSI = T - N: B as groupa, A as groupb gives A - B.
    labels = {}
    for side, prefix in (("a", "A"), ("b", "B")):
        for index, replicate in enumerate(preflight["replicates"][side], 1):
            labels[replicate] = "{}{:02d}".format(prefix, index)
    env = dict(os.environ)
    if perl5lib:
        env["PERL5LIB"] = os.pathsep.join(filter(None, [str(Path(perl5lib).resolve()), env.get("PERL5LIB")]))
    name = "PSIsigma"
    with tempfile.TemporaryDirectory(prefix="psisigma-", dir=out) as directory, \
            open(out / "psisigma.log", "w") as log:
        work = Path(directory)
        for replicate, label in labels.items():
            junction_lines, intron_lines = replicate_tables(sorted(members[replicate]), caches, universe)
            (work / (label + ".SJ.out.tab")).write_text("".join(junction_lines))
            (work / (label + ".IR.out.tab")).write_text("".join(intron_lines))
        (work / (name + ".groupa.txt")).write_text("".join(labels[r] + "\n" for r in preflight["replicates"]["b"]))
        (work / (name + ".groupb.txt")).write_text("".join(labels[r] + "\n" for r in preflight["replicates"]["a"]))
        gtf = work / "annotation.psisigma.gtf"
        shutil.copyfile(reference_dir / "annotation.psisigma.gtf", gtf)
        mapping = reference_dir / "annotation.psisigma.gtf.mapping.txt"
        database = build_database(psisigma, gtf, gtf_chromosomes(gtf), name + ".groupa.txt", name + ".groupb.txt",
                                  name, work, threads, log, env)
        found = ir_positions(database)
        if found != universe:
            missing, extra = len(set(universe) - set(found)), len(set(found) - set(universe))
            raise SystemExit("PSI-Sigma's retained introns differ from the per-run cache universe "
                             "({} missing, {} extra); irmode 0 should make them identical".format(missing, extra))
        run(["perl", Path(psisigma) / "PSIsigma-PSI-v.1.1.pl", database.name, name, settings["nread"],
             settings["skipratio"], 3, 1, 1, 0, 1, settings["irrange"], settings["variance"],
             name + ".groupa.txt", name + ".groupb.txt", settings["irclean"]], log, cwd=work, env=env)
        unfiltered = work / "{}_r{}_ir3.txt".format(name, settings["nread"])
        if not unfiltered.exists():
            raise SystemExit("PSI-Sigma's PSI step wrote no {} (see {})".format(unfiltered.name, out / "psisigma.log"))
        run(["perl", Path(psisigma) / "PSIsigma-filter-v.1.0.pl", database.name, mapping.resolve(),
             unfiltered.name, 3], log, cwd=work, env=env)
        sorted_table = work / "{}_r{}_ir3.sorted.txt".format(name, settings["nread"])
        log.flush()
        text = (out / "psisigma.log").read_text(errors="replace")
        for marker in ("(ERROR)", "Can't find", "Aborting"):
            if marker in text:
                raise SystemExit("PSI-Sigma reported '{}' (see {})".format(marker, out / "psisigma.log"))
        if not sorted_table.exists():
            raise SystemExit("PSI-Sigma's filter step wrote no {}".format(sorted_table.name))
        with tarfile.open(out / "native.tar.gz", "w:gz") as archive:
            for path in (unfiltered, sorted_table, database, work / (name + ".groupa.txt"),
                         work / (name + ".groupb.txt")):
                archive.add(str(path), arcname=path.name)
        strands = {}
        with open(mapping) as stream:
            for line in stream:
                transcript, _, strand = line.rstrip("\n").split("\t")
                strands[transcript] = strand
        rows = normalize(sorted_table, preflight, labels, base, strands)
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tested = sum(row["status"] == "tested" for row in rows)
    status, reasons = "ok", []
    if not rows:
        status, reasons = "empty", ["PSI-Sigma reported no events"]
    elif not tested:
        status, reasons = "no_tests", ["no event has a t-test p-value"]
    write_status(out / "status.json", "psisigma", status, reasons, events=len(rows), tested=tested,
                 groups={"groupa (N)": "B", "groupb (T)": "A"}, labels=labels, settings=settings,
                 statistic=STATISTIC, **base)


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def normalize(path, preflight, labels, base, strands=None):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        if len(header) != len(FILTER_HEADER):
            raise SystemExit("PSI-Sigma output has {} columns, expected {}".format(len(header), len(FILTER_HEADER)))
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(FILTER_HEADER):
                continue
            row = dict(zip(FILTER_HEADER, fields))
            database_id = row["Database ID"]
            chrom = database_id.split("_", 1)[0]
            target = row["Target Exon"]
            start = end = ""
            match = re.search(r"(\d+)-(\d+)$", target)
            if match:
                start, end = int(match.group(1)) - 1, int(match.group(2))
            dpsi, p, q = _number(row["dPSI (%)"]), _number(row["T-test p-value"]), _number(row["FDR (BH)"])
            if p is None:
                # PSI-Sigma prints dPSI 0 when it cannot test (equal means, or every replicate
                # of both groups identical): a placeholder, not an effect
                dpsi = None
            rows.append(dict(
                base, tool="psisigma", feature_id=database_id, gene_id="", gene_name=row["Gene Symbol"],
                chrom=chrom, start=start, end=end,
                coordinates=json.dumps({"event_region": row["Event Region"], "target": target,
                                        "transcript": row["Reference Transcript"], "exon_type": row["Exon Type"],
                                        "values_b": row["N Values"], "values_a": row["T Values"]}, sort_keys=True),
                strand=(strands or {}).get(re.sub(r"^(Ex|TSS)\.", "", row["Reference Transcript"]), ""),
                event_class=row["Event Type"],
                effect="" if dpsi is None else "{:.6g}".format(dpsi / 100.0), effect_definition=EFFECT,
                n_a=row["T"], n_b=row["N"], coverage="", native_statistic=row["dPSI (%)"],
                statistic_type=STATISTIC, p_value="" if p is None else "{:.6g}".format(p),
                q_value="" if q is None else "{:.6g}".format(q),
                status="tested" if p is not None and dpsi is not None else "untested"))
    return rows


def _empty_outputs(out):
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("# inference not supported\n")
    with tarfile.open(out / "native.tar.gz", "w:gz"):
        pass


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    ref = commands.add_parser("reference")
    ref.add_argument("--gtf", required=True)
    ref.add_argument("--genome", required=True)
    ref.add_argument("--psisigma", required=True, help="folder with the pinned PSI-Sigma scripts")
    ref.add_argument("--out-dir", required=True)
    ref.add_argument("--threads", type=int, default=1)
    pr = commands.add_parser("prep")
    pr.add_argument("--bam", required=True)
    pr.add_argument("--universe", required=True)
    pr.add_argument("--psisigma", required=True)
    pr.add_argument("--genome", required=True)
    pr.add_argument("--fai", default=None)
    pr.add_argument("--out-dir", required=True)
    pr.add_argument("--run-id", required=True)
    pr.add_argument("--layout", required=True, choices=("PE", "SE"))
    pr.add_argument("--threads", type=int, default=2)
    pr.add_argument("--all-junctions", action="store_true", help="keep non-canonical junctions")
    pr.add_argument("--identity", default="{}")
    cmp_ = commands.add_parser("compare")
    cmp_.add_argument("--preflight", required=True)
    cmp_.add_argument("--cache", action="append", default=[], help="RUN=DIR")
    cmp_.add_argument("--reference", required=True)
    cmp_.add_argument("--psisigma", required=True)
    cmp_.add_argument("--perl5lib", default=None)
    cmp_.add_argument("--out-dir", required=True)
    cmp_.add_argument("--threads", type=int, default=1)
    cmp_.add_argument("--ids", default="{}")
    cmp_.add_argument("--settings", default="{}")
    args = parser.parse_args(argv)
    if args.command == "reference":
        reference(args.gtf, args.genome, args.psisigma, args.out_dir, args.threads)
    elif args.command == "prep":
        prep(args.bam, args.universe, args.psisigma, args.out_dir, args.run_id, args.layout, args.genome,
             args.fai, args.threads, not args.all_junctions, json.loads(args.identity))
    else:
        compare(args.preflight, dict(item.split("=", 1) for item in args.cache), args.reference, args.psisigma,
                args.out_dir, args.threads, json.loads(args.ids), json.loads(args.settings), args.perl5lib)


if __name__ == "__main__":
    main()
