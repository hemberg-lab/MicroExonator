"""Install pinned PSI-Sigma v2.3 and the Perl modules its statistics need.

PSI-Sigma (Kuan-Ting Lin, CSHL) is free for non-commercial use by individuals
at academic or non-profit institutions; its LICENSE is kept next to the
scripts. The source is fetched at the pinned commit and every script the
workflow runs is checked against its SHA-256.

Its PSI step uses PDL, PDL::Stats, PDL::GSL::CDF and Statistics::Multtest,
which conda-forge and bioconda do not carry. They are built once with cpanm
from pinned CPAN release tarballs (SHA-256 checked) into <root>/perl5, against
the Perl and GSL of envs/umbrella-psisigma.yaml. cpanm fetches their own
dependencies from CPAN. Run this on a node with internet access.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path

COMMIT = "bfcb1423f1ec61ac48b6f92250505ec5855f4a53"   # tag v2.3
SOURCE = "https://github.com/wososa/PSI-Sigma/archive/{}.tar.gz".format(COMMIT)
SCRIPTS = {
    "PSIsigma-db-v.1.0.pl": "5d96afc7ea58f5328567c20ce8f73bcf01d5c7bc8b1461c088f41d76fcfe8654",
    "PSIsigma-ir-v.1.2.pl": "3692613c417c92d889e29a766a3476d8d8c54482ba543582ca9e00bbba296c18",
    "PSIsigma-PSI-v.1.1.pl": "286701fe7a5a2549f09f30846b196c3e2147de6c17ca2481c553fafff9b8acc3",
    "PSIsigma-filter-v.1.0.pl": "7c7218d716a8cb247a799e0d46c0ec8acd5672172a65d628da66f90f06420b1e",
    "sjFromSAMcollapseUandM_inclOverlaps.awk": "7f6d3f742cb94eac4cd1943f6c6684358c29ec9ea3eea46153eb30005599db61",
    "dummyai.pl": "432d18767dd5a6beaf1c4e31e7c1f0743bfc872878d7cefd5f777b3d9d27c610",
}
# install order: PDL-GSL and PDL-Stats build against PDL
CPAN = (
    ("PDL-2.106.tar.gz", "https://cpan.metacpan.org/authors/id/E/ET/ETJ/PDL-2.106.tar.gz",
     "f500e740d327e221985ac35f38f5e1e42bf2e3a690c0ed927ac9ed17d9d8241b"),
    ("PDL-GSL-2.103.tar.gz", "https://cpan.metacpan.org/authors/id/E/ET/ETJ/PDL-GSL-2.103.tar.gz",
     "a2e18a838589cb2e9796a04f849f9672da57b62441b720c31ad640bd16845285"),
    ("PDL-Stats-0.856.tar.gz", "https://cpan.metacpan.org/authors/id/E/ET/ETJ/PDL-Stats-0.856.tar.gz",
     "b68d3f4f8086ed357b2057d73a6592972a7d6672787a7b99499819904379540a"),
    ("Statistics-Multtest-0.14.tar.gz",
     "https://cpan.metacpan.org/authors/id/J/JO/JOKERGOO/Statistics-Multtest-0.14.tar.gz",
     "a8a33f1dbd5115071c0631c909af0dc671618baf1c5c2b63bac8e29beedc674b"),
)
MODULES = ("PDL::LiteF", "PDL::Stats", "PDL::GSL::CDF", "Statistics::Multtest")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url, target, expected=None):
    with urllib.request.urlopen(url) as response, open(target, "wb") as stream:
        shutil.copyfileobj(response, stream)
    if expected and sha256(target) != expected:
        raise SystemExit("{} does not match its pinned SHA-256".format(url))


def install_scripts(root, archive=None):
    """Unpack the pinned source into root/ and check every script."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as directory:
        if archive is None:
            archive = Path(directory) / "psisigma.tar.gz"
            fetch(SOURCE, archive)
        with tarfile.open(archive) as source:
            members = [m for m in source.getmembers() if m.isfile()]
            for member in members:
                name = member.name.split("/", 1)[-1]
                if name in SCRIPTS or name == "LICENSE":
                    with source.extractfile(member) as stream, open(root / name, "wb") as out:
                        shutil.copyfileobj(stream, out)
    bad = [name for name, digest in SCRIPTS.items() if not (root / name).exists() or sha256(root / name) != digest]
    if bad:
        raise SystemExit("PSI-Sigma scripts missing or not the pinned v2.3: {}".format(", ".join(bad)))


def install_perl(root, cpanm="cpanm", tarballs=None):
    """Build the pinned CPAN distributions into root/perl5 and check they load."""
    root = Path(root)
    lib = root / "perl5"
    with tempfile.TemporaryDirectory() as directory:
        paths = []
        for name, url, digest in CPAN:
            local = Path(tarballs) / name if tarballs else Path(directory) / name
            if not local.exists():
                fetch(url, local, digest)
            elif sha256(local) != digest:
                raise SystemExit("{} does not match its pinned SHA-256".format(local))
            paths.append(str(local))
        subprocess.run([cpanm, "--notest", "-l", str(lib)] + paths, check=True)
    env = dict(os.environ, PERL5LIB=str(lib / "lib" / "perl5"))
    check = subprocess.run(["perl"] + ["-M" + m for m in MODULES] + [
        "-e", "print join(qq{\\t}, $PDL::VERSION, $PDL::Stats::VERSION, $Statistics::Multtest::VERSION)"],
        env=env, capture_output=True, text=True)
    if check.returncode:
        raise SystemExit("PSI-Sigma's Perl modules do not load: " + check.stderr)
    versions = dict(zip(("PDL", "PDL::Stats", "Statistics::Multtest"), check.stdout.split("\t")))
    return str(lib / "lib" / "perl5"), versions


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", required=True, help="software folder (scripts, LICENSE, perl5/)")
    parser.add_argument("--record", required=True)
    parser.add_argument("--cpanm", default="cpanm")
    parser.add_argument("--source-archive", default=None, help="local copy of the pinned PSI-Sigma tarball")
    parser.add_argument("--cpan-tarballs", default=None, help="folder with the pinned CPAN tarballs")
    args = parser.parse_args(argv)
    install_scripts(args.root, args.source_archive)
    perl5lib, versions = install_perl(args.root, args.cpanm, args.cpan_tarballs)
    Path(args.record).write_text(json.dumps({
        "psisigma": "wososa/PSI-Sigma@" + COMMIT + " (v2.3)", "scripts": SCRIPTS, "perl5lib": perl5lib,
        "cpan": {name: digest for name, _, digest in CPAN}, "versions": versions,
        "licence": "free for non-commercial use by individuals at academic or non-profit institutions"},
        sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    main()
