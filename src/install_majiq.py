"""Install licensed MAJIQ v3 from user-supplied source into a workflow-owned conda env.

MAJIQ is distributed under the BioCiphers licence after registration; the
user clones or downloads majiq_academic and points `majiq_source` at it (the
folder, or a .tar.gz of it). The repository holds separate packages; as its
own Dockerfile does, this installs ./moccasin, then ./majiq (VOILA, the
viewer, is not needed by the workflow). The environment comes from
envs/umbrella-majiq.yaml at a private prefix, with HTSlib from the same
environment, and the versions are recorded. It never downloads MAJIQ.
"""

import argparse
import json
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path

PACKAGES = ("moccasin", "majiq")   # install order: majiq depends on rna_moccasin


def install(source, env_file, prefix, record, conda="conda"):
    source, prefix = Path(source), Path(prefix)
    if not source.exists():
        raise SystemExit("majiq_source not found: {}. Download MAJIQ v3 under its licence "
                         "and set majiq_source to the source directory or archive".format(source))
    if not prefix.exists():
        subprocess.run([conda, "env", "create", "-p", str(prefix), "-f", str(env_file)], check=True)
    # CMAKE_PREFIX_PATH: take compression libraries from the env, not a host without headers
    env = dict(os.environ, HTSLIB_LIBRARY_DIR=str(prefix / "lib"),
               HTSLIB_INCLUDE_DIR=str(prefix / "include"), CMAKE_PREFIX_PATH=str(prefix))
    pip = [str(prefix / "bin" / "python"), "-m", "pip", "install", "--no-cache-dir"]
    with tempfile.TemporaryDirectory() as directory:
        tree = source
        if source.is_file():
            with tarfile.open(source) as archive:
                archive.extractall(directory)
            entries = [p for p in Path(directory).iterdir() if p.is_dir()]
            tree = entries[0] if len(entries) == 1 else Path(directory)
        missing = [name for name in PACKAGES if not (tree / name / "pyproject.toml").exists()]
        if missing:
            raise SystemExit("majiq_source {} is not a majiq_academic checkout (no {}); clone "
                             "https://bitbucket.org/biociphers/majiq_academic".format(
                                 source, ", ".join(name + "/pyproject.toml" for name in missing)))
        for name in PACKAGES:
            subprocess.run(pip + [str(tree / name)], check=True, env=env)
    version = subprocess.run([str(prefix / "bin" / "majiq"), "--version"], capture_output=True,
                             text=True, env=env)
    Path(record).write_text(json.dumps({
        "source": str(source), "prefix": str(prefix), "packages": list(PACKAGES),
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
