import argparse
import gzip
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src import run_comparison_tools


class SuppaWrapperTests(unittest.TestCase):
    def test_diffsplice_gets_an_absolute_prefix(self):
        # SUPPA 2.3 merges .dpsi.temp.* from the current directory unless the
        # prefix is absolute, and then writes no .dpsi (seen on the hg38 pilot)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            joined = root / "tx.tsv.gz"
            with gzip.open(joined, "wt") as stream:
                stream.write("transcript_id\ta1\tb1\ntx1\t1\t2\n")
            preflight = {"inference_supported": True, "replicates": {"a": ["A1"], "b": ["B1"]},
                         "collapse": {"a1": "A1", "b1": "B1"}}
            args = argparse.Namespace(work=str(root / "work"), joined=str(joined), log=str(root / "log"),
                                      events=["events_SE_strict.ioe"], outputs=["suppa2/SE.dpsi"])
            commands = []
            cwd = os.getcwd()
            try:
                os.chdir(root)
                with patch.object(run_comparison_tools, "run", side_effect=lambda command, log: commands.append(command)):
                    run_comparison_tools.suppa(args, preflight)
            finally:
                os.chdir(cwd)
            diff = [command for command in commands if command[1] == "diffSplice"][0]
            prefix = diff[diff.index("-o") + 1]
            self.assertTrue(os.path.isabs(prefix), prefix)
            self.assertTrue(prefix.endswith("suppa2/SE"))


if __name__ == "__main__":
    unittest.main()
