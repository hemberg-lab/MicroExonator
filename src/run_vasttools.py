"""vast-tools module: per-run `align` caches and per-comparison differential splicing.

install-db  Once: unpack a downloaded VASTDB library, check its layout and record
            its checksum (installed.json); the archive is then deleted.
align     Per run, from the staged reads: `vast-tools align` (reads split into
          50-nt pieces, mapped to the genome and to VASTDB's junction libraries;
          strandedness detected by vast-tools). Keeps only the run's
          `to_combine` files, gzipped (vast-tools reads them gzipped), in
          to_combine.tar, plus cache.json. The split and unmapped reads are
          deleted.
compare   Per comparison, from the caches only, in a private folder:
            1. technical runs of one biological replicate are pooled with
               `vast-tools merge` (vast-tools' own way to deepen a sample);
            2. `vast-tools combine` builds the INCLUSION table;
            3. `vast-tools tidy` keeps events with enough coverage in at least
               min_samples samples of EACH group (--groups), dropping VLOW
               samples, unbalanced IR and B3 exons (--noVLOW --p_IR --noB3);
            4. `vast-tools diff` (Bayesian posterior over dPSI; MV = the |dPSI|
               exceeded with probability `prob`).
          The README recommends filtering diff with tidy, and for more than five
          samples per group, tidy + a Mann-Whitney test + |dPSI| >= 0.15
          instead (diff becomes too stringent). The call follows that:
            diff          n <= 5 per side: tidy pass and MV >= min_dpsi
            mann_whitney  n > 5 per side: tidy pass, BH q < mwu_q and
                          |mean dPSI| >= mwu_min_dpsi
          Writes
            inclusion.tsv.gz  the combine table (PSI in %, quality columns)
            tidy.tsv.gz       events passing tidy
            diff.tsv.gz       vast-tools diff (A minus B), or a not-run marker for MWU
            results.tsv.gz    normalized rows (effect = E[dPSI], A minus B)
            microexons.tsv.gz every microexon of the umbrella ME_DB with its
                              installed catalogue match, separately distinguishing
                              unreported, unquantified and absent events. None of
                              those states is contradictory differential evidence.
            status.json
"""

import argparse
import csv
import gzip
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

if __package__:
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from umbrella_modules import NORMALIZED, write_status

EFFECT = "E[dPSI] from vast-tools diff, PSI(A) - PSI(B), PSI as a fraction"
TABLES = ("inclusion.tsv.gz", "tidy.tsv.gz", "diff.tsv.gz", "results.tsv.gz", "microexons.tsv.gz")
# to_combine files vast-tools align writes per sample (IR version 2)
REQUIRED = (".eej2", ".exskX", ".micX", ".MULTI3X", ".info")
IR_FILES = {1: (".IR", ".IR.summary.txt"), 2: (".IR2", ".IR.summary_v2.txt")}
NAME = re.compile(r"[^A-Za-z0-9_]")
# vast-tools' internal species keys (bin/RunDBS_1.pl, get_internal_sp_key)
SPECIES_KEY = {"hg38": "Hs2", "hg19": "Hsa", "mm10": "Mm2", "mm9": "Mmu", "dm6": "Dme"}


def sample_name(unit_id):
    """vast-tools names samples by the file stem: no '-' or '.' (its README)."""
    name = NAME.sub("_", unit_id)
    return name if name[:1].isalpha() else "s_" + name


def vast(vast_tools, *args, cwd=None, log=None):
    command = [str(Path(vast_tools).resolve() / "vast-tools")] + [str(a) for a in args]
    if log is not None:
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
    try:
        subprocess.run(command, check=True, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        if log is None:
            raise
        log.flush()
        from collections import deque
        with open(log.name) as stream:
            tail = ''.join(deque(stream, maxlen=25))
        raise RuntimeError('vast-tools {} failed (exit {}); log: {}\n{}'.format(
            args[0], error.returncode, log.name, tail)) from error


# ------------------------------------------------------------ VASTDB

def install_db(archive, folder, species, url, record):
    """Unpack the library into `folder` (which then holds e.g. Hs2/), check it
    and write `record`; the archive is removed once the record exists."""
    import hashlib
    key = SPECIES_KEY.get(species)
    if key is None:
        raise SystemExit("no VASTDB species key for {}".format(species))
    archive, folder = Path(archive), Path(folder)
    digest = hashlib.sha256()
    with open(archive, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 22), b""):
            digest.update(block)
    with tarfile.open(archive) as library:
        names = library.getnames()
        unsafe = [n for n in names if n.startswith("/") or ".." in Path(n).parts]
        if unsafe:
            raise SystemExit("VASTDB archive has unsafe paths, e.g. {}".format(unsafe[0]))
        if not any(n.split("/")[0] == key for n in names):
            raise SystemExit("VASTDB archive {} has no {}/ folder".format(archive, key))
        library.extractall(str(folder))
    for needed in ("FILES", "TEMPLATES"):
        if not (folder / key / needed).is_dir():
            raise SystemExit("VASTDB {} has no {}/{}".format(folder, key, needed))
    version = folder / key / "VASTDB_VERSION"
    Path(record).write_text(json.dumps({
        "species": species, "species_key": key, "url": url, "archive_sha256": digest.hexdigest(),
        "archive_bytes": archive.stat().st_size,
        "vastdb_version": version.read_text().strip() if version.exists() else ""},
        sort_keys=True, indent=2) + "\n")
    archive.unlink()


# ------------------------------------------------------------ per run

def align(vast_tools, dbdir, species, reads, run_id, out_dir, threads=1, ir_version=2,
          expr=False, tmp_root=None, extra=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = sample_name(run_id)
    tmp_root = tmp_root or os.environ.get("TMPDIR") or str(out_dir)
    Path(tmp_root).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="vasttools-align-", dir=tmp_root) as directory:
        work = Path(directory)
        (work / "tmp").mkdir()
        with open(out_dir / "align.log", "w") as log:
            args = ["align"] + [str(Path(r).resolve()) for r in reads] + [
                "-sp", species, "--dbDir", str(Path(dbdir).resolve()), "-c", threads,
                "-o", str(work / "vast_out"), "--name", name, "--IR_version", ir_version,
                "--tmpDir", str(work / "tmp")]
            if expr:
                args.append("--expr")
            vast(vast_tools, *args, cwd=work, log=log)
        produced = {}
        for folder in ("to_combine", "expr_out"):
            for path in sorted((work / "vast_out" / folder).glob(name + ".*")):
                produced[folder + "/" + path.name] = path
        expected = list(REQUIRED) + list(IR_FILES[ir_version])
        missing = [suffix for suffix in expected
                   if not any(key in ("to_combine/" + name + suffix, "to_combine/" + name + suffix + ".gz")
                              for key in produced)]
        if missing:
            raise RuntimeError("vast-tools align left no {} for {} (see {})".format(
                ", ".join(missing), run_id, out_dir / "align.log"))
        eej = produced.get("to_combine/" + name + ".eej2") or produced.get("to_combine/" + name + ".eej2.gz")
        if eej.stat().st_size == 0:
            raise RuntimeError("vast-tools align wrote an empty {} for {}".format(eej.name, run_id))
        info = produced["to_combine/" + name + ".info"].read_text()
        temporary = out_dir / "to_combine.tar.tmp"
        members = {}
        with tarfile.open(temporary, "w") as archive:
            for key, path in sorted(produced.items()):
                if not path.name.endswith(".gz"):
                    with open(path, "rb") as source, gzip.open(str(path) + ".gz", "wb", compresslevel=6) as sink:
                        shutil.copyfileobj(source, sink)
                    path, key = Path(str(path) + ".gz"), key + ".gz"
                archive.add(str(path), arcname=key)
                members[key] = path.stat().st_size
    os.replace(temporary, out_dir / "to_combine.tar")
    (out_dir / "cache.json").write_text(json.dumps(dict(
        extra or {}, run_id=run_id, tool="vasttools", sample=name, species=species,
        ir_version=ir_version, expr=bool(expr), info=info.strip(), members=members),
        sort_keys=True, indent=2) + "\n")


# ------------------------------------------------------------ per comparison

def _float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def event_type(event_id, complexity):
    """exon, intron, alt_5ss or alt_3ss. vast-tools' COMPLEX column decides
    (S/C1/C2/C3/ANN/MIC*: exon; IR; Alt5; Alt3); some libraries keep raw IDs
    (e.g. FBgn..._CASSETTE1), so the ID is only a fallback."""
    complexity = complexity or ""
    if complexity.startswith("MIC") or complexity in ("S", "C1", "C2", "C3", "ANN"):
        return "exon"
    if complexity.startswith("IR"):
        return "intron"
    if complexity.startswith("Alt5"):
        return "alt_5ss"
    if complexity.startswith("Alt3"):
        return "alt_3ss"
    for key, label in (("ALTD", "alt_5ss"), ("ALTA", "alt_3ss"), ("INT", "intron"), ("EX", "exon")):
        if key in event_id:
            return label
    return "other"


def strand_of(event_id, full_coordinates, complexity=""):
    """Strand from vast-tools' FullCO column (README, combine output format):
    EX chr:C1donor,Aexon,C2acceptor (+ if C1donor < C2acceptor); ALTD
    chr:Aexon,C2acceptor; ALTA chr:C1donor,Aexon; INT chr:C1exon=C2exon:strand."""
    if ":" not in full_coordinates:
        return ""
    fields = full_coordinates.split(":")
    if event_type(event_id, complexity) == "intron":
        return fields[-1] if fields[-1] in ("+", "-") else ""
    parts = fields[1].split(",")
    first = [int(x) for x in re.findall(r"\d+", parts[0])]
    last = [int(x) for x in re.findall(r"\d+", parts[-1])]
    if not first or not last or len(parts) < 2:
        return ""
    return "+" if first[0] < last[0] else "-"


def event_class(event_id, complexity):
    if (complexity or "").startswith("MIC"):
        return "microexon"
    return {"exon": "cassette_exon", "intron": "intron_retention"}.get(
        event_type(event_id, complexity), event_type(event_id, complexity))


def read_inclusion(path):
    """{event: info} and the sample columns of a vast-tools INCLUSION table."""
    events = {}
    with open(path) as stream:
        header = stream.readline().rstrip("\n").split("\t")
        samples = header[6::2]
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            gene, event, coord, length, full, complexity = fields[:6]
            psi = {s: _float(fields[6 + 2 * i]) for i, s in enumerate(samples) if 6 + 2 * i < len(fields)}
            events[event] = dict(gene=gene, coord=coord, length=length, full=full, complexity=complexity,
                                 psi=psi)
    return events, samples


def read_table(path, key="EVENT"):
    rows = {}
    if not Path(path).exists():
        return rows
    with open(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            rows[row[key]] = row
    return rows


def read_catalog(dbdir, species):
    """Installed exon catalogue, independent of read coverage and combine output.

    combine uses FULL and rewrites old IDs via New_ID. MIC may contain legacy
    microexons absent from FULL; those are catalogued but not reported, not
    evidence against a microexon and not necessarily active quantification targets.
    """
    key = SPECIES_KEY[species]
    root = Path(dbdir) / key
    paths = [root / 'TEMPLATES' / (key + '.FULL.Template.txt.gz'),
             root / 'TEMPLATES' / (key + '.MIC.Template.txt')]
    for path in paths:
        if not path.is_file():
            raise ValueError('Missing VASTDB catalogue template: {}'.format(path))
    aliases = {}
    conversion = root / 'FILES' / ('New_ID-' + key + '.txt.gz')
    if conversion.exists():
        with gzip.open(conversion, 'rt') as stream:
            for line in stream:
                new, old, *_ = line.rstrip('\n').split('\t')
                aliases[old] = new
    entries = {}
    for path in paths:
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path, 'rt') as stream:
            header = stream.readline().rstrip('\n').split('\t')
            if header[:6] != ['GENE', 'EVENT', 'COORD', 'LENGTH', 'FullCO', 'COMPLEX']:
                raise ValueError('Unrecognized VASTDB catalogue header: {}'.format(path))
            for line in stream:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 6:
                    continue
                gene, old, coordinate, length, full, complexity = fields[:6]
                if event_type(old, complexity) != 'exon':
                    continue
                match = re.fullmatch(r'(.+):(\d+)-(\d+)', coordinate)
                if not match:
                    raise ValueError('Unrecognized exon catalogue coordinate: ' + coordinate)
                event = aliases.get(old, old)
                chrom, start, end = match.groups()
                strand = strand_of(event, full, complexity)
                entry = dict(event=event, chrom=chrom, start=int(start), end=int(end), strand=strand,
                             gene=gene, complexity=complexity)
                entries.setdefault((event, chrom, start, end, strand), entry)
    return list(entries.values())


def mann_whitney(a, b, exact_limit=20000):
    """Two-sided Mann-Whitney U p-value: exact by permutation of the rank sums
    when there are at most `exact_limit` splits, else the normal approximation
    with tie and continuity corrections."""
    if len(a) < 2 or len(b) < 2:
        return None
    values = sorted(a + b)
    ranks = {}
    position = 0
    while position < len(values):
        end = position
        while end + 1 < len(values) and values[end + 1] == values[position]:
            end += 1
        ranks[values[position]] = (position + end) / 2.0 + 1
        position = end + 1
    pooled = [ranks[v] for v in a + b]
    n_a, n_b = len(a), len(b)
    observed = sum(ranks[v] for v in a)
    expected = n_a * (n_a + n_b + 1) / 2.0
    if len(set(pooled)) == 1:
        return 1.0
    if math.comb(n_a + n_b, n_a) <= exact_limit:
        extreme = total = 0
        for chosen in itertools.combinations(range(n_a + n_b), n_a):
            total += 1
            if abs(sum(pooled[i] for i in chosen) - expected) >= abs(observed - expected) - 1e-9:
                extreme += 1
        return extreme / total
    n = n_a + n_b
    ties = sum(c ** 3 - c for c in (pooled.count(r) for r in set(pooled)))
    variance = n_a * n_b / 12.0 * ((n + 1) - ties / (n * (n - 1)))
    if variance <= 0:
        return 1.0
    z = (abs(observed - expected) - 0.5) / math.sqrt(variance)
    return min(1.0, math.erfc(max(z, 0) / math.sqrt(2)))


def benjamini_hochberg(pvalues):
    indexed = sorted((p, i) for i, p in enumerate(pvalues) if p is not None)
    q = [None] * len(pvalues)
    running = 1.0
    for rank in range(len(indexed), 0, -1):
        p, i = indexed[rank - 1]
        running = min(running, p * len(indexed) / rank)
        q[i] = running
    return q


def read_me_db(path):
    """(chrom, strand, start, end) 1-based closed, as src/insert_microexons_gtf.py reads it."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from insert_microexons_gtf import read_me_db as reader
    finally:
        sys.path.pop(0)
    return reader(path)


def compare(preflight_path, caches, vast_tools, dbdir, species, me_db, out_dir, threads=1,
            ids=None, prob=0.95, min_dpsi=0.1, min_reads=10, min_samples=2, ir_version=2,
            tidy_args="--noVLOW --p_IR --noB3", mwu_min_dpsi=0.15, mwu_q=0.05, diff_max_n=5,
            implementation="catalogue-aware-v2"):
    if implementation != 'catalogue-aware-v2':
        raise ValueError('Unsupported vast-tools comparison implementation: ' + implementation)
    preflight = json.loads(Path(preflight_path).read_text())
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ids = ids or {}
    status_path = out_dir / "status.json"
    base = dict(project_id=preflight["project_id"], comparison_id=preflight["comparison_id"], **ids)
    if not preflight["inference_supported"]:
        _empty_outputs(out_dir)
        write_status(status_path, "vasttools", "unsupported", preflight["reasons"], **base)
        return
    catalogue = read_catalog(dbdir, species)
    replicates = {side: list(preflight["replicates"][side]) for side in ("a", "b")}
    members = {}
    for run, replicate in preflight["collapse"].items():
        members.setdefault(replicate, []).append(run)
    unit_name = {run: json.loads((Path(caches[run]) / "cache.json").read_text())["sample"]
                 for run in preflight["collapse"]}
    if len(set(unit_name.values())) != len(unit_name):
        raise SystemExit("two runs of {} get the same vast-tools sample name; rename one".format(
            preflight["comparison_id"]))
    column = {}
    merges = []
    for replicate in replicates["a"] + replicates["b"]:
        runs = sorted(members[replicate])
        if len(runs) == 1:
            column[replicate] = unit_name[runs[0]]
        else:
            column[replicate] = "rep_" + sample_name(replicate)
            merges.extend((unit_name[run], column[replicate]) for run in runs)
    merged = {target for _, target in merges}
    if len(set(column.values())) != len(column) or merged & set(unit_name.values()):
        raise SystemExit("vast-tools sample names collide in {}".format(preflight["comparison_id"]))
    n_side = {side: len(replicates[side]) for side in ("a", "b")}
    samples_needed = max(1, min(min_samples, n_side["a"], n_side["b"]))
    method = "diff" if min(n_side.values()) <= diff_max_n else "mann_whitney"
    dbdir = str(Path(dbdir).resolve())
    with tempfile.TemporaryDirectory(prefix="vasttools-", dir=str(out_dir)) as directory:
        work = Path(directory)
        vast_out = work / "vast_out"
        (vast_out / "to_combine").mkdir(parents=True)
        for run in sorted(preflight["collapse"]):
            with tarfile.open(Path(caches[run]) / "to_combine.tar") as archive:
                for member in archive.getmembers():
                    if not (member.name.startswith("to_combine/") and member.isfile()):
                        continue
                    # plain files: `merge` reads only uncompressed ones
                    target = vast_out / member.name[:-3] if member.name.endswith(".gz") else vast_out / member.name
                    with archive.extractfile(member) as source, open(target, "wb") as sink:
                        shutil.copyfileobj(gzip.GzipFile(fileobj=source) if member.name.endswith(".gz")
                                           else source, sink)
        with open(out_dir / "vasttools.log", "w") as log:
            if merges:
                (work / "merge_groups.txt").write_text("".join("{}\t{}\n".format(*m) for m in merges))
                vast(vast_tools, "merge", "--groups", work / "merge_groups.txt", "-o", vast_out,
                     "--IR_version", ir_version, "--sp", species, "--dbDir", dbdir, "--move_to_PARTS",
                     cwd=work, log=log)
            vast(vast_tools, "combine", "-o", vast_out, "-sp", species, "--dbDir", dbdir,
                 "--IR_version", ir_version, "-c", threads, "--no_expr", cwd=work, log=log)
            tables = sorted(vast_out.glob("INCLUSION_LEVELS_FULL-*.tab"))
            if len(tables) != 1:
                raise SystemExit("vast-tools combine wrote {} INCLUSION tables".format(len(tables)))
            inclusion = tables[0]
            groups = "".join("{}\tA\n".format(column[r]) for r in replicates["a"]) + "".join(
                "{}\tB\n".format(column[r]) for r in replicates["b"])
            (work / "tidy_groups.txt").write_text(groups)
            vast(vast_tools, "tidy", inclusion, "--groups", work / "tidy_groups.txt", "--min_N",
                 samples_needed, "--min_SD", 0, *tidy_args.split(), "-outFile", work / "tidy.tab",
                 cwd=work, log=log)
            if method == 'diff':
                vast(vast_tools, "diff", "-a", ",".join(column[r] for r in replicates["a"]),
                     "-b", ",".join(column[r] for r in replicates["b"]), "--sampleNameA", "A",
                     "--sampleNameB", "B", "-i", inclusion, "-o", vast_out, "-d", "diff",
                     "-r", prob, "-m", min_dpsi, "-e", min_reads, "-S", samples_needed, "-c", threads,
                     "--noPDF", cwd=work, log=log)
            else:
                (vast_out / 'diff.tab').write_text('# Bayesian diff not run: using tidy + Mann-Whitney\n')
        for source, target in ((inclusion, "inclusion.tsv.gz"), (work / "tidy.tab", "tidy.tsv.gz"),
                               (vast_out / "diff.tab", "diff.tsv.gz")):
            with open(source, "rb") as stream, gzip.open(out_dir / target, "wb") as sink:
                shutil.copyfileobj(stream, sink)
        events, _ = read_inclusion(inclusion)
        tidy = read_table(work / "tidy.tab")
        diff = read_table(vast_out / "diff.tab") if method == 'diff' else {}
    names = {side: [column[r] for r in replicates[side]] for side in ("a", "b")}
    records, pvalues = [], []
    for event, info in events.items():
        values = {side: [info["psi"].get(n) for n in names[side]] for side in ("a", "b")}
        tidy_row = tidy.get(event)
        tidy_values = {side: [_float(tidy_row.get(n)) for n in names[side]] if tidy_row else []
                       for side in ("a", "b")}
        tidy_values = {side: [v for v in vs if v is not None] for side, vs in tidy_values.items()}
        p = mann_whitney(tidy_values["a"], tidy_values["b"]) if (
            tidy_row and min(len(tidy_values["a"]), len(tidy_values["b"])) >= 3) else None
        diff_row = diff.get(event, {})
        expected = _float(diff_row.get("E[dPsi]"))
        mv = _float(next((v for k, v in diff_row.items() if k and k.startswith("MV")), None))
        mean_dpsi = ((sum(tidy_values["a"]) / len(tidy_values["a"]) - sum(tidy_values["b"]) / len(tidy_values["b"]))
                     / 100.0 if tidy_values["a"] and tidy_values["b"] else None)
        chrom, start, end = "", "", ""
        match = re.match(r"^(.+):(\d+)-(\d+)$", info["coord"])
        if match:
            chrom, start, end = match.group(1), str(int(match.group(2)) - 1), match.group(3)
        pvalues.append(p)
        records.append(dict(event=event, info=info, values=values, tidy=bool(tidy_row), p=p,
                            expected=expected, mv=mv, mean_dpsi=mean_dpsi, chrom=chrom, start=start,
                            end=end, strand=strand_of(event, info["full"], info["complexity"])))
    for record, q in zip(records, benjamini_hochberg(pvalues)):
        record["q"] = q
        if method == "diff":
            tested = record["tidy"] and record["mv"] is not None
            called = tested and record["mv"] >= min_dpsi
        else:
            tested = record["tidy"] and q is not None
            called = tested and q < mwu_q and record["mean_dpsi"] is not None and abs(
                record["mean_dpsi"]) >= mwu_min_dpsi
        quantified = any(v is not None for side in record['values'].values() for v in side)
        record['status'] = ('not_quantified' if not quantified else
                            'called' if called else 'tested' if tested else 'filtered')
    normalized = []
    for record in records:
        if not any(v is not None for side in record['values'].values() for v in side):
            continue  # retain these records for catalogue mapping, not expression evidence
        info = record["info"]
        effect = record["expected"] if method == "diff" else record["mean_dpsi"]
        normalized.append({
            "tool": "vasttools", "tool_version": ids.get("tool_version", ""), "feature_id": record["event"],
            "gene_id": "", "gene_name": info["gene"], "chrom": record["chrom"], "start": record["start"],
            "end": record["end"], "coordinates": "0-based half-open alternative sequence; FullCO " + info["full"],
            "strand": record["strand"], "event_class": event_class(record["event"], info["complexity"]),
            "effect": "" if effect is None else "{:.6g}".format(effect),
            "effect_definition": EFFECT if method == "diff" else
            "mean tidy PSI(A) - mean tidy PSI(B), as a fraction",
            "n_a": sum(v is not None for v in record["values"]["a"]),
            "n_b": sum(v is not None for v in record["values"]["b"]),
            "coverage": "tidy_pass" if record["tidy"] else "tidy_fail",
            "native_statistic": "" if record["mv"] is None else "{:.6g}".format(record["mv"]),
            "statistic_type": ("vast-tools diff MV: |dPSI| exceeded with probability {}".format(prob)
                               if method == 'diff' else 'Mann-Whitney p-value; BH q-value'),
            "p_value": "" if record["p"] is None else "{:.6g}".format(record["p"]),
            "q_value": "" if record["q"] is None else "{:.6g}".format(record["q"]),
            "status": record["status"],
            **{key: base.get(key, "") for key in ("reference_id", "module_reference_id", "analysis_id",
                                                  "project_id", "comparison_id")}})
    with gzip.open(out_dir / "results.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=NORMALIZED, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(normalized)
    microexon_rows = microexon_table(read_me_db(me_db), records, method, catalogue)
    with gzip.open(out_dir / "microexons.tsv.gz", "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=MICROEXON_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(microexon_rows)
    tested = sum(r["status"] in ("tested", "called") for r in records)
    called = sum(r["status"] == "called" for r in records)
    status, reasons = "ok", []
    if not normalized:
        status, reasons = "empty", ["vast-tools quantified no event in these samples"]
    elif not tested:
        status, reasons = "no_tests", ["no event passed tidy with sufficient data for " + method]
    by_status = {}
    for row in microexon_rows:
        by_status[row["vasttools_status"]] = by_status.get(row["vasttools_status"], 0) + 1
    write_status(status_path, "vasttools", status, reasons, method=method, events=len(records),
                 tested=tested, called=called, min_samples=samples_needed, merged_replicates=sorted(
                     {m[1] for m in merges}), microexons=by_status,
                 replicates={side: replicates[side] for side in ("a", "b")}, **base)


MICROEXON_COLUMNS = ("microexon_id", "vasttools_status", "vastdb_event", "vastdb_events", "gene_name",
                     "event_class", "effect", "native_statistic", "p_value", "q_value", "method")
RANK = {"called": 0, "tested": 1, "filtered": 2, "not_quantified": 3}


def microexon_table(microexons, records, method, catalogue):
    """One row per umbrella microexon: its best-matching VASTDB exon event, by
    identical 1-based coordinates and strand, preferring the microexon (MIC)
    module. Absence is checked against installed templates, not coverage."""
    by_coordinates = {}
    catalog_coordinates = {}
    for entry in catalogue:
        catalog_coordinates.setdefault((entry['chrom'], entry['start'], entry['end']), []).append(entry)
    for record in records:
        if event_type(record["event"], record["info"]["complexity"]) != "exon" or not record["chrom"]:
            continue
        key = (record["chrom"], int(record["start"]) + 1, int(record["end"]))
        by_coordinates.setdefault(key, []).append(record)
    rows = []
    for chrom, strand, start, end in sorted(microexons):
        microexon_id = "{}_{}_{}_{}".format(chrom, strand, start - 1, end)
        candidates = by_coordinates.get((chrom, start, end), [])
        matches = [r for r in candidates if r['strand'] == strand]
        if not matches:
            entries = catalog_coordinates.get((chrom, start, end), [])
            known = [e for e in entries if e['strand'] == strand]
            known.sort(key=lambda e: (not e['complexity'].startswith('MIC'), e['event']))
            ambiguous = any(not e['strand'] for e in entries) or any(not r['strand'] for r in candidates)
            status = 'not_reported' if known else 'mapping_ambiguous' if ambiguous else 'not_in_vastdb'
            rows.append(dict(microexon_id=microexon_id, vasttools_status=status,
                             vastdb_event=known[0]['event'] if known else '',
                             vastdb_events=';'.join(sorted({e['event'] for e in known})),
                             gene_name=known[0]['gene'] if known else '',
                             event_class='microexon' if known and known[0]['complexity'].startswith('MIC') else '',
                             effect="", native_statistic="",
                             p_value="", q_value="", method=method))
            continue
        matches.sort(key=lambda r: (not r["info"]["complexity"].startswith("MIC"), RANK[r["status"]], r["event"]))
        best = matches[0]
        effect = best["expected"] if method == "diff" else best["mean_dpsi"]
        rows.append(dict(
            microexon_id=microexon_id, vasttools_status=best["status"], vastdb_event=best["event"],
            vastdb_events=";".join(r["event"] for r in matches), gene_name=best["info"]["gene"],
            event_class=event_class(best["event"], best["info"]["complexity"]),
            effect="" if effect is None else "{:.6g}".format(effect),
            native_statistic="" if best["mv"] is None else "{:.6g}".format(best["mv"]),
            p_value="" if best["p"] is None else "{:.6g}".format(best["p"]),
            q_value="" if best["q"] is None else "{:.6g}".format(best["q"]), method=method))
    return rows


def _empty_outputs(out_dir):
    for name in TABLES:
        with gzip.open(Path(out_dir) / name, "wt") as stream:
            stream.write("# inference not supported\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    db = commands.add_parser("install-db")
    db.add_argument("--archive", required=True)
    db.add_argument("--folder", required=True)
    db.add_argument("--species", required=True)
    db.add_argument("--url", default="")
    db.add_argument("--record", required=True)
    al = commands.add_parser("align")
    al.add_argument("--vast-tools", required=True, help="folder holding the vast-tools script")
    al.add_argument("--dbdir", required=True, help="VASTDB folder (holding e.g. Hs2/)")
    al.add_argument("--species", required=True, help="vast-tools assembly, e.g. hg38")
    al.add_argument("--reads", nargs="+", required=True)
    al.add_argument("--run-id", required=True)
    al.add_argument("--out-dir", required=True)
    al.add_argument("--threads", type=int, default=1)
    al.add_argument("--ir-version", type=int, default=2, choices=(1, 2))
    al.add_argument("--expr", action="store_true")
    al.add_argument("--tmp-root", default=None)
    al.add_argument("--identity", default="{}", help="JSON recorded in cache.json")
    cmp_ = commands.add_parser("compare")
    cmp_.add_argument("--preflight", required=True)
    cmp_.add_argument("--cache", action="append", default=[], help="RUN=DIR")
    cmp_.add_argument("--vast-tools", required=True)
    cmp_.add_argument("--dbdir", required=True)
    cmp_.add_argument("--species", required=True)
    cmp_.add_argument("--me-db", required=True)
    cmp_.add_argument("--out-dir", required=True)
    cmp_.add_argument("--threads", type=int, default=1)
    cmp_.add_argument("--ids", default="{}", help="JSON of identity fields for the results")
    cmp_.add_argument("--settings", default="{}", help="JSON of compare settings")
    args = parser.parse_args(argv)
    if args.command == "install-db":
        install_db(args.archive, args.folder, args.species, args.url, args.record)
    elif args.command == "align":
        align(args.vast_tools, args.dbdir, args.species, args.reads, args.run_id, args.out_dir,
              args.threads, args.ir_version, args.expr, args.tmp_root, json.loads(args.identity))
    else:
        caches = dict(item.split("=", 1) for item in args.cache)
        compare(args.preflight, caches, args.vast_tools, args.dbdir, args.species, args.me_db,
                args.out_dir, args.threads, json.loads(args.ids), **json.loads(args.settings))


if __name__ == "__main__":
    main()
