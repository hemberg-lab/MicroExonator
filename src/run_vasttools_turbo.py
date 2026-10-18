"""vast-tools-turbo: one run's `vast-tools align` split into independent chunk jobs.

split   Per run. Decides the run's strandedness exactly as `vast-tools align`
        does (first 500K reads against VASTDB's mRNA library, 0.35/0.65
        thresholds, R1 and R2 separately), then deals the read records
        round-robin into N chunk FASTQs (mates always in the same chunk).
        Writes <out>/chunks/c<k>_R<mate>.fastq.gz and <out>/split.json (the
        decision, probe fractions and read counts).
chunk   Per chunk. `vast-tools align` on one chunk with the run's strandedness
        forced (--ns, or --rc1/--nrc1 and --rc2/--nrc2), so a chunk can never
        decide differently from its siblings.
merge   Per run. `vast-tools merge` pools the chunks (vast-tools' own way to
        deepen a sample: raw and corrected counts are summed, PSI and
        complexity recomputed) and writes to_combine.tar and cache.json in the
        layout of src/run_vasttools.py `align`, so `compare` reads either.

Expected difference from a native align: the raw counts are the same; the
corrected counts (and PSI) can differ in the last rounding digits, because
vast-tools rounds per file before the sum. `merge` also overwrites the exskX complexity
column ("S" or "C", an annotation that does not depend on the reads) with "S"; that column
is restored from the chunks. Intron retention version 1 is not
supported (`merge` does not rebuild its summary).
"""

import argparse
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

if __package__:
    from src.run_vasttools import IR_FILES, REQUIRED, SPECIES_KEY, sample_name, vast
else:
    from run_vasttools import IR_FILES, REQUIRED, SPECIES_KEY, sample_name, vast

PROBE_LINES = 2000000       # RunDBS_1.pl: head -n $N for fastq input
MIN_PROBE_READS = 500       # ... "Too few first reads (<500)"
NOT_SPECIFIC = (0.35, 0.65)  # fraction of reads on the - strand of mRNA
MIN_CHUNK_READS = 1000


def chunk_name(run_id, k):
    return "{}_c{}".format(sample_name(run_id), k)


# ------------------------------------------------------------ strandedness

def probe_fractions(reads, dbdir, species, threads=1):
    """Fraction of reads mapping to the - strand of the mRNA library, per mate.

    Same command as RunDBS_1.pl: first 2M lines, bowtie -m 1 -v 2, strand column."""
    key = SPECIES_KEY[species]
    index = Path(dbdir).resolve() / key / "EXPRESSION" / "mRNA"
    fractions = []
    for path in reads:
        command = ("{unzip} {src} | head -n {n} | bowtie -q -5 0 -p {p} -m 1 -v 2 {index} - | cut -f 2".format(
            unzip="pigz -dc" if shutil.which("pigz") else "gzip -dc", src=shlex.quote(str(path)),
            n=PROBE_LINES, p=int(threads), index=shlex.quote(str(index))))
        plus = minus = 0
        process = subprocess.Popen(["bash", "-c", command], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                   text=True)
        for line in process.stdout:
            if line.rstrip("\n") == "-":
                minus += 1
            else:
                plus += 1
        process.stdout.close()
        process.wait()
        if plus + minus < MIN_PROBE_READS:
            raise RuntimeError("only {} reads of {} mapped to the mRNA library (need {}); wrong species, "
                               "or bowtie or VASTDB missing?".format(plus + minus, path, MIN_PROBE_READS))
        fractions.append(float("{:.4f}".format(minus / (plus + minus))))
    return fractions


def strand_decision(fractions):
    """RunDBS_1.pl: not strand-specific only if every mate is inside [0.35, 0.65];
    otherwise strand-specific, and a mate is reverse-complemented when its - fraction >= 0.65."""
    low, high = NOT_SPECIFIC
    if all(low <= f <= high for f in fractions):
        return {"strand_specific": False, "rc": []}
    return {"strand_specific": True, "rc": [f >= high for f in fractions]}


def strand_flags(decision):
    if not decision["strand_specific"]:
        return ["--ns"]
    flags = []
    for mate, rc in enumerate(decision["rc"], 1):
        flags.append("--{}{}".format("rc" if rc else "nrc", mate))
    return flags


# ------------------------------------------------------------ split

def _open_reader(path):
    command = ["pigz", "-dc", str(path)] if shutil.which("pigz") else ["gzip", "-dc", str(path)]
    return subprocess.Popen(command, stdout=subprocess.PIPE, bufsize=1 << 20)


def _open_writers(paths):
    zipper = ["pigz", "-1", "-p", "2"] if shutil.which("pigz") else ["gzip", "-1"]
    writers = []
    for path in paths:
        sink = open(path, "wb")
        writers.append((subprocess.Popen(zipper, stdin=subprocess.PIPE, stdout=sink, bufsize=1 << 20), sink))
    return writers


def split_reads(reads, chunks, out_dir):
    """Deal the records of `reads` (one file per mate) round-robin into `chunks` files per mate.
    Returns the reads per chunk."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mates = len(reads)
    writers = [_open_writers([out_dir / "c{}_R{}.fastq.gz".format(k, m + 1) for k in range(chunks)])
               for m in range(mates)]
    readers = [_open_reader(path) for path in reads]
    counts = [0] * chunks
    total = 0
    try:
        while True:
            records = []
            for reader in readers:
                record = [reader.stdout.readline() for _ in range(4)]
                records.append(record)
            ended = [not record[0] for record in records]
            if all(ended):
                break
            if any(ended):
                raise RuntimeError("the mates have different numbers of reads (after {} reads)".format(total))
            for record in records:
                if not record[3]:
                    raise RuntimeError("truncated FASTQ record after {} reads".format(total))
            k = total % chunks
            for m, record in enumerate(records):
                writers[m][k][0].stdin.write(b"".join(record))
            counts[k] += 1
            total += 1
    finally:
        for reader in readers:
            reader.stdout.close()
            reader.wait()
        failed = []
        for group in writers:
            for process, sink in group:
                process.stdin.close()
                if process.wait() != 0:
                    failed.append(sink.name)
                sink.close()
        if failed and sys.exc_info()[0] is None:
            raise RuntimeError("compressing chunk(s) failed: " + ", ".join(failed))
    return counts


def split(reads, dbdir, species, chunks, out_dir, run_id, threads=1, min_chunk_reads=MIN_CHUNK_READS,
          identity=None):
    out_dir = Path(out_dir)
    if chunks < 1:
        raise ValueError("chunks must be at least 1")
    fractions = probe_fractions(reads, dbdir, species, threads)
    decision = strand_decision(fractions)
    counts = split_reads(reads, chunks, out_dir / "chunks")
    total = sum(counts)
    if min(counts) < min_chunk_reads:
        raise RuntimeError("run {} has only {} reads: fewer than {} per chunk for {} chunks; use fewer "
                           "chunks (vasttools_turbo_chunks) or the native module".format(
                               run_id, total, min_chunk_reads, chunks))
    record = dict(identity or {}, run_id=run_id, mates=len(reads), chunks=chunks, reads=total,
                  reads_per_chunk=counts, probe_minus_fraction=fractions, **decision)
    (out_dir / "split.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    return record


# ------------------------------------------------------------ chunk

def count_records(path):
    reader = _open_reader(path)
    lines = 0
    for block in iter(lambda: reader.stdout.read(1 << 22), b""):
        lines += block.count(b"\n")
    reader.stdout.close()
    reader.wait()
    return lines // 4


def chunk(vast_tools, dbdir, species, split_json, index, reads, run_id, out_tar, out_meta, threads=1,
          ir_version=2, expr=False, tmp_root=None):
    """`vast-tools align` on one chunk; keeps its plain to_combine (and expr_out) files in a tar."""
    state = json.loads(Path(split_json).read_text())
    if len(reads) != state["mates"]:
        raise RuntimeError("chunk {} was given {} read files; the run has {} mates".format(
            index, len(reads), state["mates"]))
    expected = state["reads_per_chunk"][index]
    found = {count_records(path) for path in reads}
    if found != {expected}:
        raise RuntimeError("chunk {} of {} holds {} reads per mate; the split recorded {}".format(
            index, run_id, sorted(found), expected))
    name = chunk_name(run_id, index)
    out_tar, out_meta = Path(out_tar), Path(out_meta)
    out_tar.parent.mkdir(parents=True, exist_ok=True)
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out_tar.parent)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="vasttools-chunk-", dir=tmp_root) as directory:
        work = Path(directory)
        (work / "tmp").mkdir()
        log_path = work / "align.log"
        with open(log_path, "w") as log:
            args = ["align"] + [str(Path(r).resolve()) for r in reads] + [
                "-sp", species, "--dbDir", str(Path(dbdir).resolve()), "-c", threads,
                "-o", str(work / "vast_out"), "--name", name, "--IR_version", ir_version,
                "--tmpDir", str(work / "tmp")] + strand_flags(state)
            if expr:
                args.append("--expr")
            vast(vast_tools, *args, cwd=work, log=log)
        produced = {}
        for folder in ("to_combine", "expr_out"):
            for path in sorted((work / "vast_out" / folder).glob(name + ".*")):
                produced[folder + "/" + path.name] = path
        wanted = list(REQUIRED) + list(IR_FILES[ir_version])
        missing = [s for s in wanted if "to_combine/" + name + s not in produced]
        if missing:
            raise RuntimeError("vast-tools align left no {} for chunk {} of {}".format(
                ", ".join(missing), index, run_id))
        info = produced["to_combine/" + name + ".info"].read_text()
        temporary = out_tar.with_name(out_tar.name + ".tmp")
        with tarfile.open(temporary, "w") as archive:
            for key, path in sorted(produced.items()):
                archive.add(str(path), arcname=key)
    os.replace(temporary, out_tar)
    out_meta.write_text(json.dumps(dict(chunk=index, name=name, reads=expected, info=info.strip(),
                                        strand_specific=_info_is_stranded(info)),
                                   sort_keys=True, indent=2) + "\n")


def _info_is_stranded(info):
    """MergeOutputs.pl: strand-specific when the second-to-last tab field is -SS."""
    fields = info.rstrip("\n").split("\t")
    return len(fields) >= 2 and fields[-2] == "-SS"


# ------------------------------------------------------------ merge

EXSK_EVENT, EXSK_COMPLEXITY = 3, 18


def restore_exsk_complexity(merged, parts):
    """MergeOutputs.pl writes the literal "S" in the exskX complexity column, losing the
    annotation vast-tools gave each event ("S" or "C"). That value does not depend on the
    reads, so every chunk carries it unchanged; put it back and refuse if chunks disagree.
    Returns the number of events restored."""
    complexity = {}
    for part in parts:
        with open(part) as stream:
            next(stream, None)
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= EXSK_COMPLEXITY:
                    continue
                event, value = fields[EXSK_EVENT], fields[EXSK_COMPLEXITY]
                if complexity.setdefault(event, value) != value:
                    raise RuntimeError("chunks disagree on the complexity of {}: {} and {}".format(
                        event, complexity[event], value))
    restored = 0
    lines = []
    with open(merged) as stream:
        for number, line in enumerate(stream):
            fields = line.rstrip("\n").split("\t")
            if number and len(fields) > EXSK_COMPLEXITY and fields[EXSK_EVENT] in complexity \
                    and fields[EXSK_COMPLEXITY] != complexity[fields[EXSK_EVENT]]:
                fields[EXSK_COMPLEXITY] = complexity[fields[EXSK_EVENT]]
                line = "\t".join(fields) + "\n"
                restored += 1
            lines.append(line)
    with open(merged, "w") as stream:
        stream.writelines(lines)
    return restored


def merge(vast_tools, dbdir, species, split_json, tars, metas, run_id, out_dir, ir_version=2, expr=False,
          tmp_root=None, identity=None):
    if ir_version != 2:
        raise SystemExit("vasttools_turbo supports intron retention version 2 only")
    state = json.loads(Path(split_json).read_text())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = sample_name(run_id)
    if len(tars) != state["chunks"] or len(metas) != state["chunks"]:
        raise RuntimeError("{} chunks were split but {} outputs arrived".format(
            state["chunks"], min(len(tars), len(metas))))
    chunk_meta = sorted((json.loads(Path(m).read_text()) for m in metas), key=lambda c: c["chunk"])
    if [c["chunk"] for c in chunk_meta] != list(range(state["chunks"])):
        raise RuntimeError("chunk outputs do not cover chunks 0..{}".format(state["chunks"] - 1))
    if sum(state["reads_per_chunk"]) != state["reads"] or [c["reads"] for c in chunk_meta] != \
            state["reads_per_chunk"]:
        raise RuntimeError("chunk read counts do not add up to the {} reads of the run".format(state["reads"]))
    for c in chunk_meta:
        if c["strand_specific"] != state["strand_specific"]:
            raise RuntimeError("chunk {} was aligned {}strand-specific, but the run is {}".format(
                c["chunk"], "" if c["strand_specific"] else "not ", "" if state["strand_specific"] else "not "))
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out_dir)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="vasttools-merge-", dir=tmp_root) as directory:
        work = Path(directory)
        vast_out = work / "vast_out"
        for folder in ("to_combine", "expr_out"):
            (vast_out / folder).mkdir(parents=True)
        for tar_path in tars:
            with tarfile.open(tar_path) as archive:
                for member in archive.getmembers():
                    if member.isfile() and "/" in member.name and ".." not in Path(member.name).parts \
                            and member.name.split("/")[0] in ("to_combine", "expr_out"):
                        archive.extract(member, str(vast_out))
        (work / "groups.txt").write_text("".join(
            "{}\t{}\n".format(c["name"], name) for c in chunk_meta))
        with open(out_dir / "merge.log", "w") as log:
            args = ["merge", "--groups", work / "groups.txt", "-o", vast_out, "--IR_version", ir_version,
                    "--sp", species, "--dbDir", str(Path(dbdir).resolve()), "--move_to_PARTS"]
            if expr:
                args.append("--expr")
            vast(vast_tools, *args, cwd=work, log=log)
        produced = {}
        for folder in ("to_combine", "expr_out"):
            for path in sorted((vast_out / folder).glob(name + ".*")):
                produced[folder + "/" + path.name] = path
        exsk = produced.get("to_combine/" + name + ".exskX")
        if exsk is not None:
            restore_exsk_complexity(exsk, [vast_out / "to_combine" / "PARTS" / (c["name"] + ".exskX")
                                           for c in chunk_meta])
        missing = [s for s in list(REQUIRED) + list(IR_FILES[ir_version]) if "to_combine/" + name + s not in produced]
        if missing:
            raise RuntimeError("vast-tools merge left no {} for {} (see {})".format(
                ", ".join(missing), run_id, out_dir / "merge.log"))
        if produced["to_combine/" + name + ".eej2"].stat().st_size == 0:
            raise RuntimeError("vast-tools merge wrote an empty {}.eej2 for {}".format(name, run_id))
        info = produced["to_combine/" + name + ".info"].read_text()
        if _info_is_stranded(info) != state["strand_specific"]:
            raise RuntimeError("merged .info disagrees with the run's strandedness for {}".format(run_id))
        temporary = out_dir / "to_combine.tar.tmp"
        members = {}
        with tarfile.open(temporary, "w") as archive:
            for key, path in sorted(produced.items()):
                gz = Path(str(path) + ".gz")
                with open(path, "rb") as source, gzip.open(gz, "wb", compresslevel=6) as sink:
                    shutil.copyfileobj(source, sink)
                archive.add(str(gz), arcname=key + ".gz")
                members[key + ".gz"] = gz.stat().st_size
    os.replace(temporary, out_dir / "to_combine.tar")
    (out_dir / "cache.json").write_text(json.dumps(dict(
        identity or {}, run_id=run_id, tool="vasttools_turbo", sample=name, species=species,
        ir_version=ir_version, expr=bool(expr), info=info.strip(), members=members,
        chunks=state["chunks"], reads=state["reads"], reads_per_chunk=state["reads_per_chunk"],
        strand_specific=state["strand_specific"], probe_minus_fraction=state["probe_minus_fraction"]),
        sort_keys=True, indent=2) + "\n")


# ------------------------------------------------------------ CLI

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    def common(p, tool=True):
        if tool:
            p.add_argument("--vast-tools", required=True, help="folder holding the vast-tools script")
        p.add_argument("--dbdir", required=True)
        p.add_argument("--species", required=True)
        p.add_argument("--run-id", required=True)

    sp = commands.add_parser("split")
    common(sp, tool=False)
    sp.add_argument("--reads", nargs="+", required=True)
    sp.add_argument("--chunks", type=int, required=True)
    sp.add_argument("--out-dir", required=True)
    sp.add_argument("--threads", type=int, default=1)
    sp.add_argument("--min-chunk-reads", type=int, default=MIN_CHUNK_READS)
    ch = commands.add_parser("chunk")
    common(ch)
    ch.add_argument("--split-json", required=True)
    ch.add_argument("--index", type=int, required=True)
    ch.add_argument("--reads", nargs="+", required=True)
    ch.add_argument("--out-tar", required=True)
    ch.add_argument("--out-meta", required=True)
    ch.add_argument("--threads", type=int, default=1)
    ch.add_argument("--ir-version", type=int, default=2, choices=(2,))
    ch.add_argument("--expr", action="store_true")
    ch.add_argument("--tmp-root", default=None)
    mg = commands.add_parser("merge")
    common(mg)
    mg.add_argument("--split-json", required=True)
    mg.add_argument("--tars", nargs="+", required=True)
    mg.add_argument("--metas", nargs="+", required=True)
    mg.add_argument("--out-dir", required=True)
    mg.add_argument("--ir-version", type=int, default=2, choices=(2,))
    mg.add_argument("--expr", action="store_true")
    mg.add_argument("--tmp-root", default=None)
    mg.add_argument("--identity", default="{}", help="JSON recorded in cache.json")
    args = parser.parse_args(argv)
    if args.command == "split":
        split(args.reads, args.dbdir, args.species, args.chunks, args.out_dir, args.run_id, args.threads,
              args.min_chunk_reads)
    elif args.command == "chunk":
        chunk(args.vast_tools, args.dbdir, args.species, args.split_json, args.index, args.reads,
              args.run_id, args.out_tar, args.out_meta, args.threads, args.ir_version, args.expr,
              args.tmp_root)
    else:
        merge(args.vast_tools, args.dbdir, args.species, args.split_json, args.tars, args.metas,
              args.run_id, args.out_dir, args.ir_version, args.expr, args.tmp_root,
              json.loads(args.identity))


if __name__ == "__main__":
    main()
