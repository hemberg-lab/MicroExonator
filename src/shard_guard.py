"""Content identities and immutable writes for umbrella references and shards."""

import hashlib
import gzip
import json
import os
from pathlib import Path


def hisat2_members(prefix, large=False):
    extension = "ht2l" if large else "ht2"
    return ["{}.{}.{}".format(prefix, number, extension) for number in range(1, 9)]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_fixed_microexons(gtf_path, me_db_path):
    """Require every pinned 0-based ME interval as a 1-based GTF exon."""
    desired = set()
    with open(me_db_path) as stream:
        for line in stream:
            identifier = line.strip().split("\t", 1)[0]
            if not identifier or identifier.startswith("#"):
                continue
            try:
                chrom, strand, start, end = identifier.rsplit("_", 3)
                desired.add((chrom, strand, int(start), int(end)))
            except (ValueError, TypeError):
                raise ValueError("invalid fixed microexon ID: {}".format(identifier))
    if not desired:
        raise ValueError("fixed microexon universe is empty")
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    found = set()
    with opener(gtf_path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8 and fields[2] == "exon":
                coordinate = (fields[0], fields[6], int(fields[3]) - 1,
                              int(fields[4]))
                if coordinate in desired:
                    found.add(coordinate)
    missing = desired - found
    if missing:
        examples = sorted(missing)[:3]
        raise ValueError("fixed microexon universe missing {} exons from Whippet GTF; examples: {}".format(
            len(missing), examples))
    return len(desired)


def checksum_path(path, cache=None):
    """sha256 of a file, or {relative path: sha256} for a directory.

    `cache` (a dict) maps "abs path|size|mtime_ns" to a digest, so an
    unchanged file is not hashed again.
    """
    path = Path(path)
    if path.is_file():
        return _cached_sha256(path, cache)
    if not path.is_dir():
        raise ValueError("missing reference input: {}".format(path))
    members = {}
    for member in sorted(path.rglob("*")):
        if member.is_file() and member.name != ".snakemake_timestamp":
            members[str(member.relative_to(path))] = _cached_sha256(member, cache)
    if not members:
        raise ValueError("empty reference directory: {}".format(path))
    return members


def _cached_sha256(path, cache):
    if cache is None:
        return sha256_file(path)
    stat = path.stat()
    key = "{}|{}|{}".format(path.resolve(), stat.st_size, stat.st_mtime_ns)
    if key not in cache:
        cache[key] = sha256_file(path)
    return cache[key]


def reference_manifest(paths, expected_id=None, versions=None, identity=None, settings=None,
                       cache=None):
    """Hash every logical input and index member; derive the identity.

    `identity` names the entries of `paths` that define the reference_id (all
    of them when None). Built files (derived GTFs, tx2gene, indexes built by
    the workflow) are deterministic given those inputs, `versions` and
    `settings`, so they are checksummed for provenance but left out of the ID.
    That makes the ID computable before anything is built
    (`python3 src/shard_guard.py reference-id`).
    """
    checksums = {}
    for name, value in sorted(paths.items()):
        if isinstance(value, (list, tuple)):
            checksums[name] = [checksum_path(member, cache) for member in value]
        else:
            checksums[name] = checksum_path(value, cache)
    keys = sorted(paths) if identity is None else sorted(identity)
    reference_id = _identity_hash({key: checksums[key] for key in keys}, versions, settings)
    if expected_id is not None and expected_id != reference_id:
        raise ValueError("reference_id mismatch: configured {} but computed {}".format(
            expected_id, reference_id))
    return {"reference_id": reference_id, "inputs": paths, "identity_inputs": keys,
            "checksums": checksums, "versions": versions or {}, "settings": settings or {}}


def _identity_hash(checksums, versions, settings):
    identity = {"checksums": checksums, "versions": versions or {}}
    if settings:
        identity["settings"] = settings
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()[:16]


REFERENCE_FILES = {"genome": "genome_fasta", "annotation": "annotation_gtf",
                   "whippet_annotation": "whippet_gtf", "me_db": "me_db",
                   "salmon_gtf": "salmon_gtf", "transcripts": "transcriptome_fasta",
                   "decoys": "decoys", "splice_sites": "splice_sites",
                   "annotation_bed12": "annotation_bed12"}


def reference_identity(reference, insert_microexons=True):
    """(paths, settings) that define a reference_id, from config `umbrella_reference`.

    Configured files and any prebuilt index count; built outputs do not.
    """
    paths = {name: reference[key] for name, key in REFERENCE_FILES.items() if reference.get(key)}
    large = reference.get("hisat2_index_type", "small") == "large"
    if reference.get("hisat2_index_prefix"):
        paths["hisat2"] = hisat2_members(reference["hisat2_index_prefix"], large)
    if reference.get("whippet_index"):
        paths["whippet"] = [reference["whippet_index"], reference["whippet_index"] + ".exons.tab.gz"]
    if reference.get("salmon_index"):
        paths["salmon"] = reference["salmon_index"]
    settings = {"hisat2_index_type": "large" if large else "small",
                "insert_microexons": bool(insert_microexons)}
    return paths, settings


def config_reference_id(reference, insert_microexons=True, cache_path=None):
    """The reference_id of config `umbrella_reference`, from configured inputs only.

    With `cache_path` (JSON), file digests are reused while a file's path,
    size and mtime are unchanged: only the first call hashes the genome.
    """
    cache = {}
    if cache_path and Path(cache_path).is_file():
        try:
            cache = json.loads(Path(cache_path).read_text())
        except ValueError:
            cache = {}
    before = dict(cache)
    paths, settings = reference_identity(reference, insert_microexons)
    reference_id = reference_manifest(paths, versions=reference.get("versions", {}),
                                      settings=settings, cache=cache)["reference_id"]
    if cache_path and cache != before:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        temporary = Path("{}.{}.tmp".format(cache_path, os.getpid()))
        temporary.write_text(json.dumps(cache, sort_keys=True, indent=1) + "\n")
        os.replace(temporary, cache_path)   # atomic: parallel parses never see half a file
    return reference_id


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def write_immutable(path, content):
    path = Path(path)
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError("immutable output differs: {}".format(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(str(path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
    except FileExistsError:
        if path.read_bytes() != content:
            raise ValueError("immutable output differs: {}".format(path))
        return
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(content)


def write_immutable_bundle(outputs, guard_path, metadata):
    """Preflight all outputs and their identity before making any new writes."""
    guard_path = Path(guard_path)
    checksums = {str(path): hashlib.sha256(data).hexdigest()
                 for path, data in sorted(outputs.items())}
    guard = _json_bytes({"metadata": metadata, "sha256": checksums})
    for path, data in outputs.items():
        path = Path(path)
        if path.exists() and path.read_bytes() != data:
            raise ValueError("immutable output differs: {}".format(path))
    if guard_path.exists() and guard_path.read_bytes() != guard:
        raise ValueError("immutable guard differs: {}".format(guard_path))
    for path, data in outputs.items():
        write_immutable(path, data)
    write_immutable(guard_path, guard)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description="Print the umbrella reference_id of a config.")
    parser.add_argument("command", choices=["reference-id"])
    parser.add_argument("--configfile", required=True, help="Snakemake config (YAML or JSON)")
    args = parser.parse_args(argv)
    with open(args.configfile) as stream:
        text = stream.read()
    try:
        config = json.loads(text)
    except ValueError:
        import yaml
        config = yaml.safe_load(text)
    insert = str(config.get("umbrella_insert_microexons", True)).lower() not in ("false", "f", "0", "no")
    print(config_reference_id(config["umbrella_reference"], insert))


if __name__ == "__main__":
    main()
