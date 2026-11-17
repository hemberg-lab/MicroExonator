"""Core --circ quantification and a checked import into the BS-analysis cache."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess

if __package__:
    from src.whippet_flags import quant_flags
    from src.run_whippet_circ import count_lines, whippet_version
else:
    from whippet_flags import quant_flags
    from run_whippet_circ import count_lines, whippet_version

KEPT = ("psi.gz", "jnc.gz", "gene.tpm.gz", "isoform.tpm.gz")


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def core_quant(reads, index, prefix, run_id, reference_id, julia="julia", whippet_bin="", flags="", processing_sha256=None):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    flags = shlex.split(quant_flags(flags))
    command = [julia, str(Path(whippet_bin) / "whippet-quant.jl"), *reads,
               "-x", str(index), "-o", str(prefix), *flags]
    subprocess.run(command, check=True)
    for suffix in (*KEPT, "map.gz"):
        if not Path(str(prefix) + "." + suffix).is_file():
            raise ValueError("Whippet wrote no quant." + suffix)
    events, circular = count_lines(str(prefix) + ".psi.gz")
    record = dict(schema="whippet-core-v1", run_id=run_id, reference_id=reference_id,
                  whippet_version=whippet_version(whippet_bin), flags=flags, circ=True, processing_sha256=processing_sha256,
                  events=events, bs_events=circular,
                  sha256={suffix: digest(str(prefix) + "." + suffix) for suffix in KEPT})
    Path(str(prefix) + ".provenance.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    return record


def validate_core(source, run_id, identity, verify_hashes=True):
    source = Path(source)
    marker = source / "quant.provenance.json"
    if not marker.is_file():
        raise ValueError("Missing core --circ provenance; backfill legacy Whippet first: " + str(marker))
    record = json.loads(marker.read_text())
    if record.get("schema") != "whippet-core-v1" or not record.get("circ") or "--circ" not in record.get("flags", []):
        raise ValueError("Core provenance does not certify --circ: " + str(marker))
    if record.get("run_id") != run_id:
        raise ValueError("Core provenance belongs to another run: " + str(marker))
    if record.get("reference_id") != identity.get("reference_id"):
        raise ValueError("Core reference does not match requested cache")
    if "flags" in identity and record.get("flags") != identity["flags"]:
        raise ValueError("Core flags do not match requested cache settings")
    if "processing_sha256" in identity and record.get("processing_sha256") != identity["processing_sha256"]:
        raise ValueError("Core processing identity does not match requested cache")
    for suffix in KEPT:
        path = source / ("quant." + suffix)
        if not path.is_file():
            raise ValueError("Missing retained core output: " + str(path))
        if verify_hashes and digest(path) != record.get("sha256", {}).get(suffix):
            raise ValueError("Core checksum mismatch: " + str(path))
    return record


def import_core(source, out_dir, run_id, identity):
    record = validate_core(source, run_id, identity)
    out = Path(out_dir)
    # Never replace a cache that might belong to a previous computation.
    if out.exists() and any(out.iterdir()):
        raise ValueError("Cache destination is not empty: " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    for suffix in KEPT:
        origin, target = Path(source) / ("quant." + suffix), out / ("quant." + suffix)
        try:
            os.link(origin, target)  # no extra copy of the retained quantifications
        except OSError:
            shutil.copy2(origin, target)
    result = dict(identity, run_id=run_id, tool="whippet_circ", source="core --circ",
                  whippet_version=record["whippet_version"], events=record["events"],
                  bs_events=record["bs_events"], core_provenance=record)
    (out / "cache.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    core = subs.add_parser("quant")
    core.add_argument("--reads", nargs="+", required=True)
    core.add_argument("--index", required=True)
    core.add_argument("--prefix", required=True)
    core.add_argument("--reference-id", required=True)
    core.add_argument("--processing-sha256", required=True)
    core.add_argument("--julia", default="julia")
    core.add_argument("--whippet-bin", required=True)
    core.add_argument("--flags", default="")
    reuse = subs.add_parser("import")
    reuse.add_argument("--source", required=True)
    reuse.add_argument("--out-dir", required=True)
    reuse.add_argument("--identity", required=True)
    for sub in (core, reuse):
        sub.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if args.command == "quant":
        core_quant(args.reads, args.index, args.prefix, args.run_id, args.reference_id,
                   args.julia, args.whippet_bin, args.flags, args.processing_sha256)
    else:
        import_core(args.source, args.out_dir, args.run_id, json.loads(args.identity))


if __name__ == "__main__":
    main()
