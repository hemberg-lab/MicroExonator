"""Check the NumPy-dependent delta wrapper's supported and marker paths."""

import argparse
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))
import run_comparison_tools  # noqa: E402


class MicroexonatorDeltaWrapperTests(unittest.TestCase):
    def test_unsupported_comparison_writes_marker_without_running_tool(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "delta.tsv"
            args = argparse.Namespace(outputs=[str(output)], work=str(Path(tmp) / "work"))
            with patch.object(run_comparison_tools.subprocess, "run") as run:
                run_comparison_tools.microexonator(args, {"inference_supported": False})
            self.assertEqual(output.read_text(), "# inference not supported\n")
            run.assert_not_called()

    def test_supported_comparison_uses_conda_python_and_cleans_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            output = directory / "delta.tsv"
            work = directory / "work"
            args = argparse.Namespace(outputs=[str(output)], work=str(work),
                                      joined="joined.tsv.gz", microexons="microexons.tsv",
                                      options="--min-reads 10")
            preflight = {"inference_supported": True,
                         "runs": {"a": ["run_a"], "b": ["run_b"]}}
            with patch.object(run_comparison_tools, "write_delta_inputs") as rebuild, \
                 patch.object(run_comparison_tools.subprocess, "run") as run:
                run_comparison_tools.microexonator(args, preflight)
            rebuild.assert_called_once()
            command = run.call_args.args[0]
            self.assertEqual(command[:2], ["python3", "src/me_delta.py"])
            self.assertEqual(command[-2:], ["--min-reads", "10"])
            self.assertIn("run_a.tsv.gz", command[command.index("-a") + 1])
            self.assertIn("run_b.tsv.gz", command[command.index("-b") + 1])
            self.assertFalse(work.exists())


if __name__ == "__main__":
    unittest.main()
