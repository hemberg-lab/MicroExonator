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
                    (Path(cwd) / "leafcutter_ds_effect_sizes.txt").write_text("intron\tdeltapsi_b\n")
                    (Path(cwd) / "comparison_perind_numers.counts.gz").write_bytes(b"introns")

            result = root / "leafcutter_cluster_significance.txt"
            exons = root / "exons.txt.gz"
            exons.write_bytes(b"")
            with patch("src.run_leafcutter.subprocess.run", side_effect=fake_command):
                run(preflight, joined, result, root / "leafcutter.log", threads=4,
                    effect_sizes=root / "effects.txt", introns=root / "introns.counts.gz", exons=exons)
            self.assertEqual([command[0] for command in commands],
                             ["leafcutter-cluster", "leafcutter-ds"])
            self.assertIn("--baseline_group", commands[1])
            ds = commands[1]
            self.assertEqual(ds[ds.index("--min_samples_per_intron") + 1], "1")
            self.assertEqual(ds[ds.index("--min_samples_per_group") + 1], "1")
            self.assertIn("p.adjust", result.read_text())
            self.assertEqual(ds[ds.index("--exon_file") + 1], str(exons.resolve()))
            self.assertEqual((root / "effects.txt").read_text(), "intron\tdeltapsi_b\n")
            self.assertEqual((root / "introns.counts.gz").read_bytes(), b"introns")
            self.assertFalse(list(root.glob("leafcutter-*")))

    def test_sample_limits_follow_smallest_group(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            joined = root / "junctions.tsv.gz"
            names = ["A1", "A2", "A3", "B1", "B2", "B3", "B4", "B5", "B6"]
            with gzip.open(joined, "wt") as output:
                output.write("chrom\tstart\tend\tstrand\t" + "\t".join(names) + "\n")
                output.write("chr1\t21\t30\t+\t" + "\t".join(["5"] * len(names)) + "\n")
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({"inference_supported": True,
                                             "replicates": {"a": names[:3], "b": names[3:]}}))
            commands = []

            def fake_command(command, cwd, **kwargs):
                commands.append(command)
                if command[0] == "leafcutter-ds":
                    (Path(cwd) / "leafcutter_ds_cluster_significance.txt").write_text("cluster\n")

            with patch("src.run_leafcutter.subprocess.run", side_effect=fake_command):
                run(preflight, joined, root / "out.txt", root / "leafcutter.log")
            ds = commands[1]
            self.assertEqual(ds[ds.index("--min_samples_per_intron") + 1], "3")
            self.assertEqual(ds[ds.index("--min_samples_per_group") + 1], "3")
            self.assertIn("calibrated down to 4", (root / "leafcutter.log").read_text())

    def test_unsupported_comparison_does_not_invoke_tools(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = root / "preflight.json"
            preflight.write_text('{"inference_supported": false}')
            result = root / "result.txt"
            with patch("src.run_leafcutter.subprocess.run") as command:
                run(preflight, root / "missing.tsv.gz", result, root / "leafcutter.log",
                    effect_sizes=root / "effects.txt", introns=root / "introns.counts.gz")
            command.assert_not_called()
            self.assertEqual(result.read_text(), "# inference not supported\n")
            self.assertEqual((root / "effects.txt").read_text(), "# inference not supported\n")


if __name__ == "__main__":
    unittest.main()
