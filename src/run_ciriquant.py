"""CIRIquant module: per-run circRNA counts against a frozen catalogue, and comparison tables.

index    Once per catalogue: CIRIquant's generate_index + build_index (the same pseudo-circular
         FASTA and hisat2-build call it makes per run), kept as a shared reference. prep links
         it in place of the per-run build; counts are identical (native test).
prep     Per run, while its staged FASTQs and the shared umbrella HISAT2 BAM exist:
         CIRIquant 1.1.3's own steps, in CIRIquant's order (circ.proc): pseudo-circular
         FASTA from the catalogue, HISAT2 index, realignment of all read pairs, BSJ reads
         (proc_denovo_bam), FSJ reads and BSJ re-checks on the genome BAM (proc_genome_bam),
         library sizes (bam_stat). Keeps counts.tsv.gz for every catalogue circRNA, not just
         those CIRIquant's GTF would print (it skips BSJ = 0), and cache.json.
compare  Per comparison, from the kept counts only: per-replicate BSJ, FSJ and CIRIquant's
         junction ratio 2*BSJ/(2*BSJ + FSJ); A minus B of the replicate means. Descriptive:
         no differential test is claimed for these ratios in this release.

Why not the `CIRIquant` command: version 1.1.3 still has Python 2 code on paths every run
uses, and its config check requires BWA, Java and Perl indexes that the --bam/--bed route
never touches. This module sets CIRIquant's tool and reference globals directly and swaps
only these helpers (PATCHES), leaving alignment, filtering and counting native:
  load_fai     float file offsets under Python 3 (int(...) * int(...) / int(...))
  grouper      list * float chunk size once a catalogue has more than 500 circRNAs per thread
(denovo_worker passes a float junction position to pysam get_overlap, which pysam 0.22+
accepts, so it stays native.)
format_output (sorted(..., cmp=...)) and check_samtools_version (cmp) are not called.
Counting units: BSJ and FSJ are read pairs (CIRIquant's query_prefix), as in CIRIquant.
"""

import argparse
import gzip
import json
import logging
import multiprocessing
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

PATCHES = ("load_fai: integer byte offsets", "grouper: integer chunk size")
CIRIQUANT_VERSION = "1.1.3"
UNITS = "read pairs (CIRIquant query_prefix), mates of one fragment counted once"
COUNT_COLUMNS = ("circ_id", "chrom", "start", "end", "strand", "bsj", "fsj", "junction_ratio", "measurement")


def library_type(declared, summary=None):
    """CIRIquant -l code and the strandedness it stands for.

    0 unstranded; 1 read1 on the transcript's strand (umbrella 'secondstrand');
    2 read1 on the opposite strand (umbrella 'firststrand', dUTP). 'auto' reads the
    umbrella junction summary's inferred strandedness, like the other umbrella tools."""
    if declared == "auto":
        if not summary:
            raise ValueError("strandedness auto needs the run's junction summary")
        declared = json.loads(Path(summary).read_text())["inferred_strandedness"]
        if declared == "undetermined":
            declared = "unstranded"
    codes = {"unstranded": 0, "secondstrand": 1, "firststrand": 2}
    if declared not in codes:
        raise ValueError("unknown strandedness {!r}".format(declared))
    return codes[declared], declared


def junction_ratio(bsj, fsj):
    """CIRIquant's junction ratio, or None when there is no support at all."""
    denominator = 2.0 * bsj + fsj
    return 2.0 * bsj / denominator if denominator > 0 else None


def _patch(circ):
    from itertools import zip_longest

    def load_fai(fname):
        faidx = {}
        with open(fname, 'r') as f:
            for line in f:
                content = line.rstrip().split('\t')
                chrom, length, start, eff_length, line_length = content
                shift_length = int(length) * int(line_length) // int(eff_length)
                faidx[chrom] = [int(start), shift_length]
        return faidx

    def grouper(iterable, n, fillvalue=None):
        args = [iter(iterable)] * int(n)
        return zip_longest(*args, fillvalue=None)

    circ.load_fai = load_fai
    circ.grouper = grouper


def _open_out(path):
    return gzip.open(path, "wt")


def _ciriquant(genome):
    """CIRIquant's modules, patched, with its tool and reference globals set."""
    try:
        multiprocessing.set_start_method("fork", force=True)
    except (RuntimeError, ValueError):
        pass
    from CIRIquant import circ, utils
    _patch(circ)
    utils.HISAT2 = shutil.which("hisat2")
    utils.SAMTOOLS = shutil.which("samtools")
    # not resolve(): the module genome may be a link whose .fai sits beside the link
    utils.FASTA = os.path.abspath(genome)
    if not utils.HISAT2 or not utils.SAMTOOLS:
        raise SystemExit("hisat2 and samtools must be on PATH")
    if not Path(utils.FASTA + ".fai").exists():
        raise SystemExit("genome {} needs its samtools .fai next to it".format(genome))
    return circ, utils


def _logger(log_file):
    logger = logging.getLogger("CIRIquant")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)


def index(catalogue, genome, out_dir, threads=1):
    """Shared pseudo-circular index: out_dir/pseudo.{fa,*.ht2}, indexed.txt, index.json."""
    circ, _ = _ciriquant(genome)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log_file = str((out / "index.log").resolve())
    _logger(log_file)
    circ_info = circ.load_bed(str(Path(catalogue).resolve()))
    fasta = out / "pseudo.fa"
    circ.generate_index(log_file, circ_info, str(fasta))
    with open(fasta) as stream:
        indexed = sorted(line[1:].strip() for line in stream if line.startswith(">"))
    # build_index writes <outdir>/circ/<prefix>_index; build into a scratch folder, then move
    with tempfile.TemporaryDirectory(prefix="ciriquant-index-", dir=out) as directory:
        (Path(directory) / "circ").mkdir()
        circ.build_index(log_file, threads, str(fasta), directory, "pseudo")
        built = sorted((Path(directory) / "circ").glob("pseudo_index.*.ht2*"))
        if not built:
            raise SystemExit("hisat2-build wrote no index (see {})".format(log_file))
        for path in built:
            os.replace(path, out / path.name)
    (out / "indexed.txt").write_text("".join(cid + "\n" for cid in indexed))
    record = {"catalogue_sha256": _sha256(catalogue), "indexed": len(indexed), "ciriquant_version": CIRIQUANT_VERSION,
              "patches": list(PATCHES), "files": sorted(path.name for path in built)}
    (out / "index.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def _sha256(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prep(reads, bam, catalogue, genome, out_dir, run_id, threads=1, anchor=5, strandedness="unstranded",
         summary=None, layout="PE", identity=None, tmp_root=None, index_dir=None):
    if layout != "PE" or len(reads) != 2:
        raise SystemExit("CIRIquant needs paired-end reads (its realignment uses -1/-2); run {} is {}".format(
            run_id, layout))
    lib_type, strand_name = library_type(strandedness, summary)
    circ, utils = _ciriquant(genome)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log_file = str((out / "ciriquant.log").resolve())
    _logger(log_file)
    reads = [str(Path(r).resolve()) for r in reads]
    bam = str(Path(bam).resolve())
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    prefix = "run"
    with tempfile.TemporaryDirectory(prefix="ciriquant-", dir=tmp_root) as directory:
        work = Path(directory)
        circ_dir = work / "circ"
        circ_dir.mkdir()
        circ_info = circ.load_bed(str(Path(catalogue).resolve()))
        if index_dir:
            # the shared index CIRIquant's own functions built once for this catalogue
            record = json.loads((Path(index_dir) / "index.json").read_text())
            if record["catalogue_sha256"] != _sha256(catalogue):
                raise SystemExit("pseudo-circular index {} was built from another catalogue".format(index_dir))
            for name in record["files"]:
                (circ_dir / name.replace("pseudo_index", prefix + "_index")).symlink_to(
                    (Path(index_dir) / name).resolve())
            indexed = set((Path(index_dir) / "indexed.txt").read_text().split())
        else:
            circ_fasta = str(circ_dir / (prefix + "_index.fa"))
            circ.generate_index(log_file, circ_info, circ_fasta)
            with open(circ_fasta) as stream:
                indexed = {line[1:].strip() for line in stream if line.startswith(">")}
            circ.build_index(log_file, threads, circ_fasta, str(work), prefix)
        denovo_bam = circ.denovo_alignment(log_file, threads, reads, str(work), prefix)
        candidates = circ.proc_denovo_bam(denovo_bam, threads, circ_info, anchor, lib_type)
        bsj_reads, fsj_reads = circ.proc_genome_bam(bam, threads, circ_info, candidates, anchor, str(circ_dir), False)
        total_reads, mapped_reads = circ.bam_stat(bam)
        native = circ.expression_level(circ_info, bsj_reads, fsj_reads)
    rows = []
    for chrom in sorted(circ_info):
        for cid, parser in sorted(circ_info[chrom].items(), key=lambda item: (item[1].start, item[1].end)):
            if cid not in indexed:
                rows.append((cid, chrom, parser.start, parser.end, parser.strand, "", "", "", "not_indexed"))
                continue
            bsj = len(bsj_reads[cid]) if cid in bsj_reads else 0
            fsj = len(fsj_reads[cid]) if cid in fsj_reads else 0
            if cid in native and (native[cid]["bsj"], native[cid]["fsj"]) != (bsj, fsj):
                raise SystemExit("count export differs from CIRIquant's expression_level for " + cid)
            ratio = junction_ratio(bsj, fsj)
            rows.append((cid, chrom, parser.start, parser.end, parser.strand, bsj, fsj,
                         "" if ratio is None else "{:.6f}".format(ratio), "measured" if ratio is not None else "no_support"))
    with _open_out(out / "counts.tsv.gz.tmp") as stream:
        stream.write("\t".join(COUNT_COLUMNS) + "\n")
        for row in rows:
            stream.write("\t".join(str(value) for value in row) + "\n")
    os.replace(out / "counts.tsv.gz.tmp", out / "counts.tsv.gz")
    circular_reads = 2 * sum(len(bsj_reads[cid]) for cid in bsj_reads)
    record = dict(identity or {}, run_id=run_id, tool="ciriquant", ciriquant_version=CIRIQUANT_VERSION,
                  patches=list(PATCHES), library_type=lib_type, strandedness=strand_name, anchor=anchor,
                  counting_units=UNITS, total_reads=total_reads, mapped_reads=mapped_reads,
                  circular_reads=circular_reads, catalogue_entries=len(rows),
                  measured=sum(1 for row in rows if row[-1] == "measured"),
                  pseudo_index="shared" if index_dir else "per_run",
                  not_indexed=sum(1 for row in rows if row[-1] == "not_indexed"))
    (out / "cache.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    return record


def read_counts(path):
    counts = {}
    with gzip.open(path, "rt") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        for line in stream:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            counts[row["circ_id"]] = row
    return counts


def compare(preflight, caches, out_dir, ids=None, settings=None, selection=None):
    """Descriptive per-comparison tables from kept counts (A minus B)."""
    from src.umbrella_modules import NORMALIZED, write_status
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    settings = dict({"min_bsj": 1}, **(settings or {}))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ids = ids or {}
    base = dict(ids, project_id=preflight.get("project_id", ""), comparison_id=preflight.get("comparison_id", ""),
                selection=selection or {})
    reasons = []
    if not all(preflight["replicates"][side] for side in ("a", "b")):
        _empty(out)
        write_status(out / "status.json", "ciriquant", "unsupported", ["a side has no replicate"], **base)
        return
    records = {}
    for unit, cache in caches.items():
        record = json.loads((Path(cache) / "cache.json").read_text())
        if record.get("run_id") != unit:
            raise SystemExit("CIRIquant cache {} belongs to {}, not {}".format(cache, record.get("run_id"), unit))
        records[unit] = record
    references = {record.get("module_reference_id") for record in records.values()}
    if len(references) > 1:
        raise SystemExit("CIRIquant caches come from different catalogues or settings: {}".format(sorted(references)))
    members = defaultdict(list)
    for unit, replicate in preflight["collapse"].items():
        members[replicate].append(unit)
    side_of = {replicate: side for side in ("a", "b") for replicate in preflight["replicates"][side]}
    replicates = [r for side in ("a", "b") for r in preflight["replicates"][side]]
    bsj, fsj, info = defaultdict(dict), defaultdict(dict), {}
    library = {}
    for replicate in replicates:
        totals = [0, 0]
        for unit in sorted(members[replicate]):
            counts = read_counts(Path(caches[unit]) / "counts.tsv.gz")
            totals[0] += records[unit]["total_reads"]
            totals[1] += records[unit]["mapped_reads"]
            for cid, row in counts.items():
                info.setdefault(cid, row)
                if row["measurement"] == "not_indexed":
                    bsj[cid].setdefault(replicate, None)
                    fsj[cid].setdefault(replicate, None)
                    continue
                if bsj[cid].get(replicate, 0) is None:
                    continue
                bsj[cid][replicate] = bsj[cid].get(replicate, 0) + int(row["bsj"])
                fsj[cid][replicate] = fsj[cid].get(replicate, 0) + int(row["fsj"])
        library[replicate] = {"side": side_of[replicate], "units": ",".join(sorted(members[replicate])),
                              "total_reads": totals[0], "mapped_reads": totals[1]}
    if len(set(preflight["collapse"].values())) != len(preflight["collapse"]):
        reasons.append("technical units summed within biological replicates")
    _matrices(out, replicates, bsj, fsj, library)
    rows = []
    for cid in sorted(bsj):
        values = {r: (bsj[cid].get(r), fsj[cid].get(r)) for r in replicates}
        total_bsj = sum(b for b, _ in values.values() if b is not None)
        if total_bsj < settings["min_bsj"]:
            continue
        ratios = {side: [junction_ratio(*values[r]) for r in preflight["replicates"][side]
                         if values[r][0] is not None and junction_ratio(*values[r]) is not None]
                  for side in ("a", "b")}
        n_a, n_b = len(ratios["a"]), len(ratios["b"])
        effect = (sum(ratios["a"]) / n_a - sum(ratios["b"]) / n_b) if n_a and n_b else ""
        row = info[cid]
        rows.append(dict(base, tool="ciriquant", tool_version="CIRIquant " + CIRIQUANT_VERSION, feature_id=cid,
                         gene_id="", gene_name="", chrom=row["chrom"], start=row["start"], end=row["end"],
                         coordinates=cid, strand=row["strand"], event_class="back_splice_junction",
                         effect="" if effect == "" else "{:.6f}".format(effect),
                         effect_definition="mean CIRIquant junction ratio 2*BSJ/(2*BSJ+FSJ) per replicate, A minus B",
                         n_a=n_a, n_b=n_b, coverage=total_bsj, native_statistic="", statistic_type="descriptive",
                         p_value="", q_value="", status="descriptive" if effect != "" else "undefined_ratio"))
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("\t".join(NORMALIZED) + "\n")
        for row in rows:
            stream.write("\t".join(str(row.get(column, "")) for column in NORMALIZED) + "\n")
    if not rows:
        write_status(out / "status.json", "ciriquant", "empty", reasons + ["no catalogue circRNA reached min_bsj"],
                     **base)
        return
    reasons.append("descriptive junction ratios; CIRIquant's own test (CIRI_DE_replicate) is in native_de/")
    write_status(out / "status.json", "ciriquant", "no_tests", reasons, circrnas=len(rows), **base)


def _matrices(out, replicates, bsj, fsj, library):
    for name, table in (("bsj_counts", bsj), ("fsj_counts", fsj)):
        with gzip.open(out / (name + ".tsv.gz"), "wt") as stream:
            stream.write("\t".join(["circ_id"] + replicates) + "\n")
            for cid in sorted(table):
                stream.write("\t".join([cid] + ["NA" if table[cid].get(r) is None else str(table[cid].get(r))
                                                for r in replicates]) + "\n")
    with open(out / "library.tsv", "w") as stream:
        stream.write("replicate\tside\tunits\ttotal_reads\tmapped_reads\n")
        for replicate in replicates:
            entry = library[replicate]
            stream.write("\t".join([replicate, entry["side"], entry["units"], str(entry["total_reads"]),
                                    str(entry["mapped_reads"])]) + "\n")


def _empty(out):
    from src.umbrella_modules import NORMALIZED
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("\t".join(NORMALIZED) + "\n")
    for name in ("bsj_counts", "fsj_counts"):
        with gzip.open(out / (name + ".tsv.gz"), "wt") as stream:
            stream.write("circ_id\n")
    (out / "library.tsv").write_text("replicate\tside\tunits\ttotal_reads\tmapped_reads\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prep")
    p.add_argument("--reads", nargs="+", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--catalogue", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--run-id", required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--anchor", type=int, default=5)
    p.add_argument("--strandedness", default="unstranded")
    p.add_argument("--summary")
    p.add_argument("--layout", default="PE")
    p.add_argument("--identity", default="{}")
    p.add_argument("--tmp-root")
    p.add_argument("--index-dir", help="shared pseudo-circular index from the index command")
    i = commands.add_parser("index")
    i.add_argument("--catalogue", required=True)
    i.add_argument("--genome", required=True)
    i.add_argument("--out-dir", required=True)
    i.add_argument("--threads", type=int, default=1)
    c = commands.add_parser("compare")
    c.add_argument("--preflight", required=True)
    c.add_argument("--cache", action="append", default=[], help="unit=cache_dir")
    c.add_argument("--out-dir", required=True)
    c.add_argument("--ids", default="{}")
    c.add_argument("--settings", default="{}")
    c.add_argument("--selection", default="{}")
    args = parser.parse_args(argv)
    if args.command == "prep":
        prep(args.reads, args.bam, args.catalogue, args.genome, args.out_dir, args.run_id, args.threads, args.anchor,
             args.strandedness, args.summary, args.layout, json.loads(args.identity), args.tmp_root, args.index_dir)
    elif args.command == "index":
        print(json.dumps(index(args.catalogue, args.genome, args.out_dir, args.threads)))
    else:
        caches = dict(item.split("=", 1) for item in args.cache)
        compare(args.preflight, caches, args.out_dir, json.loads(args.ids), json.loads(args.settings),
                json.loads(args.selection))


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    sys.exit(main())
