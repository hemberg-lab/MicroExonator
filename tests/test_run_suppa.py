import argparse
import gzip
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src import run_comparison_tools


# Stands in for SUPPA 2.3: psiPerEvent writes <o>.psi; diffSplice writes
# <o>.dpsi.temp.0 with conditions named after the .psi files, then merges
# every *.dpsi.temp.* in the prefix's folder (the current directory when the
# prefix is relative), whatever their prefix, and deletes them, as
# lib/diff_tools.merge_temp_output_files does. FAKE_SUPPA_EXTRA_COLUMN makes
# diffSplice emit a stray column.
FAKE_SUPPA = r'''
import os
import sys

args = sys.argv[1:]
out = args[args.index("-o") + 1]
if args[0] == "psiPerEvent":
    with open(out + ".psi", "w") as stream:
        stream.write("s1\n")
    sys.exit(0)
psi = args[args.index("-p") + 1:args.index("-p") + 3]
cond = "-".join(os.path.basename(p)[:-len(".psi")] for p in psi)
kind = os.path.basename(psi[0]).split("_")[0]
with open(out + ".dpsi.temp.0", "w") as stream:
    stream.write("Event_id\t{0}_dPSI\t{0}_p-val\n".format(cond))
    stream.write("g;{}:1-2\t0.1\t0.5\n".format(kind))
folder = os.path.dirname(out) if os.path.isabs(out) else os.getcwd()
temps = sorted(f for f in os.listdir(folder) if ".dpsi.temp." in f)
header, rows = [], []
for name in temps:
    with open(os.path.join(folder, name)) as stream:
        lines = stream.read().splitlines()
    header += lines[0].split("\t")[1:]
    rows += lines[1:]
if os.environ.get("FAKE_SUPPA_EXTRA_COLUMN"):
    header.append("stray")
with open(out + ".dpsi", "w") as stream:
    stream.write("\t".join(header) + "\n" + "".join(r + "\n" for r in rows))
for name in temps:
    os.remove(os.path.join(folder, name))
'''


class SuppaWrapperTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        fake = bin_dir / "suppa.py"
        fake.write_text("#!" + sys.executable + "\n" + FAKE_SUPPA)
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
        joined = self.root / "tx.tsv.gz"
        with gzip.open(joined, "wt") as stream:
            stream.write("transcript_id\ta1\tb1\ntx1\t1\t2\n")
        self.preflight = {"inference_supported": True, "replicates": {"a": ["A1"], "b": ["B1"]},
                          "collapse": {"a1": "A1", "b1": "B1"}}
        self.kinds = ("SE", "A5", "AL")
        self.args = argparse.Namespace(
            work=str(self.root / "work"), joined=str(joined), log=str(self.root / "log"),
            events=["events_{}_strict.ioe".format(k) for k in self.kinds],
            outputs=["suppa2/{}.dpsi".format(k) for k in self.kinds])
        self.env = patch.dict(os.environ, {"PATH": str(bin_dir) + os.pathsep + os.environ["PATH"]})
        self.env.start()
        self.cwd = os.getcwd()
        os.chdir(self.root)

    def tearDown(self):
        os.chdir(self.cwd)
        self.env.stop()
        self.tmp.cleanup()

    def test_each_dpsi_holds_only_its_own_event_type(self):
        # the hg38 pilot: a run with relative prefixes left every type's
        # .dpsi.temp.0 in suppa2/, and the next run's SE merge swallowed them
        (self.root / "suppa2").mkdir()
        for kind in ("AL", "RI", "MX"):
            (self.root / "suppa2" / (kind + ".dpsi.temp.0")).write_text(
                "Event_id\t{0}_b-{0}_a_dPSI\t{0}_b-{0}_a_p-val\ng;{0}:9-9\t0\t1\n".format(kind))
        (self.root / "leftover.dpsi.temp.0").write_text("Event_id\tx_dPSI\tx_p-val\n")
        run_comparison_tools.suppa(self.args, self.preflight)
        for kind in self.kinds:
            lines = (self.root / "suppa2" / (kind + ".dpsi")).read_text().splitlines()
            self.assertEqual(lines[0], "{0}_b-{0}_a_dPSI\t{0}_b-{0}_a_p-val".format(kind))
            self.assertEqual(lines[1:], ["g;{}:1-2\t0.1\t0.5".format(kind)])
        self.assertFalse((self.root / "work").exists())

    def test_diffsplice_prefix_is_private_and_absolute(self):
        commands = []
        real_run = run_comparison_tools.run
        with patch.object(run_comparison_tools, "run",
                          side_effect=lambda c, log: (commands.append(c), real_run(c, log))):
            run_comparison_tools.suppa(self.args, self.preflight)
        prefixes = [c[c.index("-o") + 1] for c in commands if c[1] == "diffSplice"]
        self.assertEqual(len(set(os.path.dirname(p) for p in prefixes)), len(self.kinds))
        for kind, prefix in zip(self.kinds, prefixes):
            self.assertTrue(os.path.isabs(prefix), prefix)
            self.assertEqual(Path(prefix).parent.name, kind)

    def test_unexpected_header_fails_loudly(self):
        with patch.dict(os.environ, {"FAKE_SUPPA_EXTRA_COLUMN": "1"}):
            with self.assertRaises(SystemExit) as caught:
                run_comparison_tools.suppa(self.args, self.preflight)
        self.assertIn("SE", str(caught.exception.code))
        self.assertFalse((self.root / "suppa2" / "SE.dpsi").exists())


if __name__ == "__main__":
    unittest.main()
