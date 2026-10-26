"""circRNA discovery: one shared STAR chimeric alignment per run, CIRCexplorer2 and DCC.

refflat        umbrella annotation GTF -> CIRCexplorer2's gene annotation (refFlat layout:
               gene, transcript, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,
               exonStarts, exonEnds; 0-based starts).
star           Per run, from the staged reads: STAR 2.7.11b with chimeric detection only
               (--outSAMtype None, --chimOutType Junctions, legacy junction format), with
               DCC's recommended joint-mapping settings, which CIRCexplorer2 also accepts.
               Gzipped reads are streamed through named pipes this script feeds, not
               STAR's --readFilesCommand (STAR spawns a script without a #! line, which
               works on Linux only). Keeps Chimeric.out.junction.gz, SJ.out.tab.gz,
               Log.final.out and cache.json.
circexplorer2  Per run, from the STAR cache only: CIRCexplorer2 parse -t STAR, annotate.
               Keeps its native tables and calls.tsv.gz (circ_id chrom:start|end, 1-based,
               the catalogue's convention).
support        Per comparison, CIRCexplorer2 calls of its runs: reads per replicate, number
               of replicates per side detecting each circRNA, catalogue support. Detection
               only, no test.
dcc            Per comparison, from its runs' STAR caches: DCC 0.5.0 detection and its
               replicate filter (-Nr), unstranded mode with the annotation GTF for strand
               and gene. DCC runs its "CircSkip" step whenever -an is given; that step is
               switched off with an empty -C file (DCC reads -C only for host-gene
               counting, -G, which is off), because it needs annotation features this
               route does not provide and is not used here.
Technical units of one biological replicate are pooled before DCC and summed in support.
"""

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

STAR_CHIMERIC = ("--chimSegmentMin", "15", "--chimJunctionOverhangMin", "15", "--chimScoreMin", "15",
                 "--chimScoreSeparation", "10", "--alignSJoverhangMin", "15", "--alignSJDBoverhangMin", "15",
                 "--outSJfilterOverhangMin", "15", "15", "15", "15", "--chimMultimapNmax", "0",
                 "--chimOutJunctionFormat", "0", "--chimOutType", "Junctions", "--outSAMtype", "None")
STAR_FILES = ("Chimeric.out.junction.gz", "SJ.out.tab.gz", "Log.final.out", "cache.json")
CALL_COLUMNS = ("circ_id", "chrom", "start", "end", "strand", "reads", "gene", "circ_type", "exon_count")


def circ_id(chrom, start1, end):
    return "{}:{}|{}".format(chrom, start1, end)


def gtf_to_refflat(gtf, out):
    exons = defaultdict(list)
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            transcript = re.search(r'transcript_id "([^"]+)"', fields[8])
            gene = re.search(r'gene_name "([^"]+)"', fields[8]) or re.search(r'gene_id "([^"]+)"', fields[8])
            if not transcript or not gene:
                continue
            exons[(gene.group(1), transcript.group(1), fields[0], fields[6])].append((int(fields[3]) - 1, int(fields[4])))
    with open(out, "w") as stream:
        for (gene, transcript, chrom, strand), blocks in sorted(exons.items(), key=lambda item: (item[0][2], min(item[1]), item[0][1])):
            blocks.sort()
            stream.write("\t".join([gene, transcript, chrom, strand, str(blocks[0][0]), str(blocks[-1][1]),
                                    str(blocks[0][0]), str(blocks[-1][1]), str(len(blocks)),
                                    "".join("{},".format(s) for s, _ in blocks),
                                    "".join("{},".format(e) for _, e in blocks)]) + "\n")
    return len(exons)


def _run(command, log, cwd=None):
    log.write("$ " + " ".join(str(part) for part in command) + "\n")
    log.flush()
    result = subprocess.run([str(part) for part in command], cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit("{} failed with exit {} (see {})".format(Path(str(command[0])).name, result.returncode,
                                                                   log.name))


def _log_final(path):
    stats = {}
    for line in Path(path).read_text().splitlines():
        if "|" in line:
            key, value = (part.strip() for part in line.split("|", 1))
            stats[key] = value
    keys = {"input_reads": "Number of input reads", "unique_percent": "Uniquely mapped reads %",
            "chimeric_reads": "Number of chimeric reads", "chimeric_percent": "% of chimeric reads"}
    return {name: stats.get(key, "") for name, key in keys.items()}


def star(reads, index, out_dir, run_id, threads=1, identity=None, tmp_root=None, star_bin="STAR"):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    version = subprocess.run([star_bin, "--version"], capture_output=True, text=True).stdout.strip()
    with tempfile.TemporaryDirectory(prefix="circrna-star-", dir=tmp_root) as directory, \
            open(out / "star.log", "w") as log:
        work = Path(directory)
        inputs, feeders = [], []
        for number, path in enumerate(reads, 1):
            path = str(Path(path).resolve())
            if path.endswith(".gz"):
                fifo = work / "R{}.fastq".format(number)
                os.mkfifo(fifo)
                feeders.append(subprocess.Popen("gzip -dc {} > {}".format(_quote(path), _quote(str(fifo))), shell=True))
                inputs.append(str(fifo))
            else:
                inputs.append(path)
        try:
            _run([star_bin, "--runThreadN", threads, "--genomeDir", str(Path(index).resolve()),
                  "--readFilesIn"] + inputs + ["--outFileNamePrefix", str(work) + "/"] + list(STAR_CHIMERIC), log)
        finally:
            for feeder in feeders:
                feeder.wait()
        if any(feeder.returncode for feeder in feeders):
            raise SystemExit("decompressing the reads of {} failed".format(run_id))
        for name in ("Chimeric.out.junction", "SJ.out.tab"):
            with open(work / name, "rb") as source, gzip.open(out / (name + ".gz.tmp"), "wb") as sink:
                shutil.copyfileobj(source, sink)
            os.replace(out / (name + ".gz.tmp"), out / (name + ".gz"))
        shutil.copyfile(work / "Log.final.out", out / "Log.final.out")
    stats = _log_final(out / "Log.final.out")
    record = dict(identity or {}, run_id=run_id, tool="circrna_star", star_version=version,
                  parameters=list(STAR_CHIMERIC), mates=len(reads), **stats)
    (out / "cache.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    return record


def _quote(text):
    import shlex
    return shlex.quote(text)


def junction_lines(path):
    """Chimeric.out.junction(.gz) data lines, without STAR's header or comment lines."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith(("#", "chr_donorA")) or not line.strip():
                continue
            yield line


def circexplorer2(star_cache, refflat, genome, out_dir, run_id, identity=None, tmp_root=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="circexplorer2-", dir=tmp_root) as directory, \
            open(out / "circexplorer2.log", "w") as log:
        work = Path(directory)
        junctions = work / "Chimeric.out.junction"
        lines = 0
        with open(junctions, "w") as sink:
            for line in junction_lines(Path(star_cache) / "Chimeric.out.junction.gz"):
                sink.write(line)
                lines += 1
        _run(["CIRCexplorer2", "parse", "-t", "STAR", junctions, "-b", work / "back_spliced_junction.bed"], log)
        _run(["CIRCexplorer2", "annotate", "-r", Path(refflat).resolve(), "-g", Path(genome).resolve(),
              "-b", work / "back_spliced_junction.bed", "-o", work / "circularRNA_known.txt"], log)
        for name in ("back_spliced_junction.bed", "circularRNA_known.txt"):
            with open(work / name, "rb") as source, gzip.open(out / (name + ".gz"), "wb") as sink:
                shutil.copyfileobj(source, sink)
        calls = read_circexplorer2(work / "circularRNA_known.txt")
    with gzip.open(out / "calls.tsv.gz", "wt") as stream:
        stream.write("\t".join(CALL_COLUMNS) + "\n")
        for row in calls:
            stream.write("\t".join(str(row[column]) for column in CALL_COLUMNS) + "\n")
    record = dict(identity or {}, run_id=run_id, tool="circexplorer2", circexplorer2_version="2.3.8",
                  chimeric_lines=lines, circrnas=len(calls), reads=sum(row["reads"] for row in calls))
    (out / "cache.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    return record


def read_circexplorer2(path):
    """CIRCexplorer2 annotate output (BED12+; column 13 readNumber) -> call rows."""
    rows = []
    with open(path) as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15:
                continue
            chrom, start1, end = fields[0], int(fields[1]) + 1, int(fields[2])
            rows.append({"circ_id": circ_id(chrom, start1, end), "chrom": chrom, "start": start1, "end": end,
                         "strand": fields[5], "reads": int(fields[12]), "gene": fields[14], "circ_type": fields[13],
                         "exon_count": fields[9]})
    rows.sort(key=lambda row: (row["chrom"], row["start"], row["end"], row["strand"]))
    return rows


def _catalogue_ids(path):
    if not path:
        return None
    with open(path) as stream:
        return {line.split("\t")[3] for line in stream if line.strip()}


def _members(preflight):
    members = defaultdict(list)
    for unit, replicate in preflight["collapse"].items():
        members[replicate].append(unit)
    return members


def _write_results(out, rows, tool, version):
    from src.umbrella_modules import NORMALIZED
    with gzip.open(out / "results.tsv.gz", "wt") as stream:
        stream.write("\t".join(NORMALIZED) + "\n")
        for row in rows:
            stream.write("\t".join(str(row.get(column, "")) for column in NORMALIZED) + "\n")


def _support_rows(preflight, per_replicate, info, catalogue, base, tool, version, definition):
    """per_replicate: {circ_id: {replicate: reads}} -> normalized rows and support rows."""
    replicates = {side: preflight["replicates"][side] for side in ("a", "b")}
    rows, support = [], []
    for cid in sorted(per_replicate):
        counts = per_replicate[cid]
        n_a = sum(1 for r in replicates["a"] if counts.get(r, 0) > 0)
        n_b = sum(1 for r in replicates["b"] if counts.get(r, 0) > 0)
        chrom, start, end, strand, gene = info[cid]
        supported = "" if catalogue is None else ("yes" if cid in catalogue else "no")
        rows.append(dict(base, tool=tool, tool_version=version, feature_id=cid, gene_id="", gene_name=gene,
                         chrom=chrom, start=start, end=end, coordinates=cid, strand=strand,
                         event_class="back_splice_junction", effect="", effect_definition=definition,
                         n_a=n_a, n_b=n_b, coverage=sum(counts.values()), native_statistic="",
                         statistic_type="detection", p_value="", q_value="", status="detected"))
        support.append([cid, chrom, start, end, strand, gene, supported, n_a, n_b] +
                       [counts.get(r, 0) for side in ("a", "b") for r in replicates[side]])
    header = (["circ_id", "chrom", "start", "end", "strand", "gene", "catalogue_supported", "replicates_a",
               "replicates_b"] + [r for side in ("a", "b") for r in replicates[side]])
    return rows, header, support


def support(preflight, caches, out_dir, ids=None, catalogue=None, selection=None):
    """CIRCexplorer2 calls of one comparison: detection support per replicate."""
    from src.umbrella_modules import write_status
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    base = dict(ids or {}, project_id=preflight.get("project_id", ""), comparison_id=preflight.get("comparison_id", ""),
                selection=selection or {})
    per_replicate, info = defaultdict(lambda: defaultdict(int)), {}
    for replicate, units in _members(preflight).items():
        for unit in units:
            record = json.loads((Path(caches[unit]) / "cache.json").read_text())
            if record.get("run_id") != unit:
                raise SystemExit("CIRCexplorer2 cache {} belongs to {}, not {}".format(caches[unit], record.get("run_id"), unit))
            with gzip.open(Path(caches[unit]) / "calls.tsv.gz", "rt") as stream:
                header = stream.readline().rstrip("\n").split("\t")
                for line in stream:
                    row = dict(zip(header, line.rstrip("\n").split("\t")))
                    per_replicate[row["circ_id"]][replicate] += int(row["reads"])
                    info.setdefault(row["circ_id"], (row["chrom"], row["start"], row["end"], row["strand"], row["gene"]))
    rows, header, table = _support_rows(preflight, per_replicate, info, _catalogue_ids(catalogue), base,
                                        "circexplorer2", "CIRCexplorer2 2.3.8",
                                        "detection only: CIRCexplorer2 BSJ reads per replicate, no effect size")
    _write_results(out, rows, "circexplorer2", "CIRCexplorer2 2.3.8")
    _write_support(out, header, table)
    if not rows:
        write_status(out / "status.json", "circexplorer2", "empty", ["no circRNA called in these runs"], **base)
    else:
        write_status(out / "status.json", "circexplorer2", "no_tests", ["detection only: no differential test"],
                     circrnas=len(rows), **base)


def _write_support(out, header, table):
    with gzip.open(out / "support.tsv.gz", "wt") as stream:
        stream.write("\t".join(header) + "\n")
        for row in table:
            stream.write("\t".join(str(value) for value in row) + "\n")


def dcc(preflight, caches, gtf, out_dir, ids=None, settings=None, catalogue=None, selection=None, threads=1,
        tmp_root=None):
    from src.umbrella_modules import write_status
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    settings = dict({"min_reads": 2, "min_replicates": 1, "min_length": 30, "max_length": 1000000, "chrM": True},
                    **(settings or {}))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    base = dict(ids or {}, project_id=preflight.get("project_id", ""), comparison_id=preflight.get("comparison_id", ""),
                selection=selection or {})
    members = _members(preflight)
    replicates = [r for side in ("a", "b") for r in preflight["replicates"][side]]
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="dcc-", dir=tmp_root) as directory, open(out / "dcc.log", "w") as log:
        work = Path(directory)
        sheet = []
        for index, replicate in enumerate(replicates):
            # DCC names samples by file name; a plain index keeps them parseable
            path = work / "s{:03d}.junction".format(index)
            with open(path, "w") as sink:
                for unit in sorted(members[replicate]):
                    record = json.loads((Path(caches[unit]) / "cache.json").read_text())
                    if record.get("run_id") != unit:
                        raise SystemExit("STAR cache {} belongs to {}, not {}".format(caches[unit], record.get("run_id"), unit))
                    for line in junction_lines(Path(caches[unit]) / "Chimeric.out.junction.gz"):
                        sink.write(line)
            sheet.append(str(path))
        (work / "samplesheet").write_text("\n".join(sheet) + "\n")
        (work / "no_skip.circ").write_text("")
        command = ["DCC", "@" + str(work / "samplesheet"), "-D", "-N", "-an", Path(gtf).resolve(), "-F",
                   "-Nr", settings["min_reads"], settings["min_replicates"], "-n", settings["min_length"],
                   "-m", settings["max_length"], "-C", work / "no_skip.circ", "-T", threads,
                   "-O", str(work / "out") + "/", "-t", str(work / "tmp") + "/"]
        if settings["chrM"]:
            command.insert(command.index("-F") + 1, "-M")
        _run(command, log, cwd=work)
        counts = (work / "out" / "CircRNACount").read_text().splitlines()
        coordinates = (work / "out" / "CircCoordinates").read_text().splitlines()
        for name, lines in (("CircRNACount", counts), ("CircCoordinates", coordinates)):
            with gzip.open(out / (name + ".tsv.gz"), "wt") as stream:
                stream.write("\n".join(lines) + "\n")
    per_replicate, info = defaultdict(dict), {}
    labels = {"s{:03d}.junction".format(i): r for i, r in enumerate(replicates)}
    header = counts[0].split("\t")
    for count_line, coordinate_line in zip(counts[1:], coordinates[1:]):
        values = count_line.split("\t")
        coordinate = coordinate_line.split("\t")
        if values[:3] != coordinate[:3]:
            raise SystemExit("DCC CircRNACount and CircCoordinates rows disagree: {} vs {}".format(values[:3], coordinate[:3]))
        chrom, start1, end = coordinate[0], int(coordinate[1]), int(coordinate[2])
        cid = circ_id(chrom, start1, end)
        info[cid] = (chrom, start1, end, coordinate[5], coordinate[3])
        for column, value in zip(header[3:], values[3:]):
            per_replicate[cid][labels[column]] = int(value)
    rows, support_header, table = _support_rows(preflight, per_replicate, info, _catalogue_ids(catalogue), base,
                                                "dcc", "DCC 0.5.0",
                                                "detection only: DCC chimeric reads per replicate after -Nr filtering")
    _write_results(out, rows, "dcc", "DCC 0.5.0")
    _write_support(out, support_header, table)
    reasons = ["detection only: no differential test (DCC/CircTest linear counts are not part of this release)"]
    if len(set(preflight["collapse"].values())) != len(preflight["collapse"]):
        reasons.append("technical units pooled within biological replicates before DCC")
    write_status(out / "status.json", "dcc", "no_tests" if rows else "empty",
                 reasons if rows else ["DCC kept no circRNA after filtering"], circrnas=len(rows), settings=settings,
                 **base)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    r = commands.add_parser("refflat")
    r.add_argument("--gtf", required=True)
    r.add_argument("--out", required=True)
    s = commands.add_parser("star")
    s.add_argument("--reads", nargs="+", required=True)
    s.add_argument("--index", required=True)
    s.add_argument("--out-dir", required=True)
    s.add_argument("--run-id", required=True)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--identity", default="{}")
    s.add_argument("--tmp-root")
    c = commands.add_parser("circexplorer2")
    c.add_argument("--star-cache", required=True)
    c.add_argument("--refflat", required=True)
    c.add_argument("--genome", required=True)
    c.add_argument("--out-dir", required=True)
    c.add_argument("--run-id", required=True)
    c.add_argument("--identity", default="{}")
    c.add_argument("--tmp-root")
    for name in ("support", "dcc"):
        p = commands.add_parser(name)
        p.add_argument("--preflight", required=True)
        p.add_argument("--cache", action="append", default=[], help="unit=cache_dir")
        p.add_argument("--out-dir", required=True)
        p.add_argument("--ids", default="{}")
        p.add_argument("--catalogue")
        p.add_argument("--selection", default="{}")
        if name == "dcc":
            p.add_argument("--gtf", required=True)
            p.add_argument("--settings", default="{}")
            p.add_argument("--threads", type=int, default=1)
            p.add_argument("--tmp-root")
    args = parser.parse_args(argv)
    if args.command == "refflat":
        print(gtf_to_refflat(args.gtf, args.out), "transcripts")
    elif args.command == "star":
        star(args.reads, args.index, args.out_dir, args.run_id, args.threads, json.loads(args.identity), args.tmp_root)
    elif args.command == "circexplorer2":
        circexplorer2(args.star_cache, args.refflat, args.genome, args.out_dir, args.run_id,
                      json.loads(args.identity), args.tmp_root)
    else:
        caches = dict(item.split("=", 1) for item in args.cache)
        if args.command == "support":
            support(args.preflight, caches, args.out_dir, json.loads(args.ids), args.catalogue, json.loads(args.selection))
        else:
            dcc(args.preflight, caches, args.gtf, args.out_dir, json.loads(args.ids), json.loads(args.settings),
                args.catalogue, json.loads(args.selection), args.threads, args.tmp_root)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    sys.exit(main())
