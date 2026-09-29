"""Exercise the Snakemake-facing LeafCutter wrapper without installing PyTorch."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src.run_leafcutter import run


class RunLeafcutterTest(unittest.TestCase):
    def test_managed_commands_use_generated_junctions_and_groups(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            joined = root / "junctions.tsv.gz"
            with gzip.open(joined, "wt") as output:
                output.write("chrom\tstart\tend\tstrand\tA\tB\n")
                output.write("chr1\t21\t30\t+\t40\t20\n")
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({"inference_supported": True,
                                             "replicates": {"a": ["A"], "b": ["B"]}}))
            commands = []

            def fake_command(command, cwd, **kwargs):
                commands.append(command)
                if command[0] == "leafcutter-cluster":
                    self.assertEqual((Path(cwd) / "A.junc").read_text(),
                                     "chr1\t19\t31\t.\t40\t+\t19\t31\t255,0,0\t2\t1,1\t0,11\n")
                    self.assertEqual((Path(cwd) / "groups.txt").read_text(),
                                     "A\ta\nB\tb\n")
                else:
                    (Path(cwd) / "leafcutter_ds_cluster_significance.txt").write_text(
                        "cluster\tp.adjust\nchr1:clu_1_+\t0.01\n")

            result = root / "leafcutter_cluster_significance.txt"
            with patch("src.run_leafcutter.subprocess.run", side_effect=fake_command):
                run(preflight, joined, result, root / "leafcutter.log", threads=4)
            self.assertEqual([command[0] for command in commands],
                             ["leafcutter-cluster", "leafcutter-ds"])
            self.assertIn("--baseline_group", commands[1])
            self.assertIn("p.adjust", result.read_text())
            self.assertFalse(list(root.glob("leafcutter-*")))

    def test_unsupported_comparison_does_not_invoke_tools(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = root / "preflight.json"
            preflight.write_text('{"inference_supported": false}')
            result = root / "result.txt"
            with patch("src.run_leafcutter.subprocess.run") as command:
                run(preflight, root / "missing.tsv.gz", result, root / "leafcutter.log")
            command.assert_not_called()
            self.assertEqual(result.read_text(), "# inference not supported\n")


if __name__ == "__main__":
    unittest.main()
