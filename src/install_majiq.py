"""Install licensed MAJIQ v3 from user-supplied source into a workflow-owned conda env.

MAJIQ is distributed under the BioCiphers licence after registration; the
user downloads the source and points `majiq_source` at it (a directory or a
.tar.gz). This creates the environment from envs/umbrella-majiq.yaml at a
private prefix, installs moccasin and MAJIQ with pip inside it (HTSlib from
the same environment), and records the versions. It never downloads MAJIQ.
"""

import argparse
import json
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path

MOCCASIN = "git+https://bitbucket.org/biociphers/moccasin@new_moccasin"


def install(source, env_file, prefix, record, conda="conda"):
    source, prefix = Path(source), Path(prefix)
    if not source.exists():
        raise SystemExit("majiq_source not found: {}. Download MAJIQ v3 under its licence "
                         "and set majiq_source to the source directory or archive".format(source))
    if not prefix.exists():
        subprocess.run([conda, "env", "create", "-p", str(prefix), "-f", str(env_file)], check=True)
    env = dict(os.environ, HTSLIB_LIBRARY_DIR=str(prefix / "lib"),
               HTSLIB_INCLUDE_DIR=str(prefix / "include"))
    pip = [str(prefix / "bin" / "python"), "-m", "pip", "install", "--no-cache-dir"]
    subprocess.run(pip + [MOCCASIN], check=True, env=env)
    with tempfile.TemporaryDirectory() as directory:
        tree = source
        if source.is_file():
            with tarfile.open(source) as archive:
                archive.extractall(directory)
            entries = [p for p in Path(directory).iterdir() if p.is_dir()]
            tree = entries[0] if len(entries) == 1 else Path(directory)
        subprocess.run(pip + [str(tree)], check=True, env=env)
    version = subprocess.run([str(prefix / "bin" / "majiq"), "--version"], capture_output=True,
                             text=True, env=env)
    Path(record).write_text(json.dumps({
        "source": str(source), "prefix": str(prefix), "moccasin": MOCCASIN,
        "majiq_version": (version.stdout or version.stderr).strip()}, sort_keys=True, indent=2) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--env-file", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--record", required=True)
    parser.add_argument("--conda", default="conda")
    args = parser.parse_args(argv)
    install(args.source, args.env_file, args.prefix, args.record, args.conda)


if __name__ == "__main__":
    main()
