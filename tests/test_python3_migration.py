import os
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]


class Python3MigrationTests(unittest.TestCase):
    def test_all_source_scripts_compile_with_python3(self):
        result = subprocess.run(
            [sys.executable, "-m", "compileall", "-q", "src"],
            cwd=REPOSITORY,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        self.assertEqual(result.returncode, 0, result.stdout)

    def test_gtf_to_bed12_runs_with_python3(self):
        gtf = (
            'chr1\tsource\ttranscript\t101\t300\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "tx1";\n'
            'chr1\tsource\texon\t101\t150\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "tx1";\n'
            'chr1\tsource\texon\t251\t300\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "tx1";\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            gtf_path = pathlib.Path(directory) / "input.gtf"
            gtf_path.write_text(gtf)
            result = subprocess.run(
                [sys.executable, "src/GTFtoBED12.py", str(gtf_path)],
                cwd=REPOSITORY,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            result.stdout,
            "chr1\t100\t300\ttx1\t0\t+\t100\t300\t0\t2\t50,50\t0,150\n",
        )

    def test_workflow_explicitly_invokes_python3(self):
        workflow_files = [
            *REPOSITORY.glob("*.smk"),
            *REPOSITORY.glob("rules/*.smk"),
            REPOSITORY / "config.py",
            REPOSITORY / "src" / "Snakefile",
        ]
        python2_files = [
            path.relative_to(REPOSITORY)
            for path in workflow_files
            if "python2" in path.read_text()
        ]
        implicit_python_files = [
            path.relative_to(REPOSITORY)
            for path in workflow_files
            if re.search(r"\bpython\s+src/", path.read_text())
        ]

        self.assertEqual(python2_files, [])
        self.assertEqual(implicit_python_files, [])

    def test_pipeline_environments_do_not_pin_python2(self):
        python2_environments = [
            path.relative_to(REPOSITORY)
            for path in (REPOSITORY / "envs").glob("*.yaml")
            if "python=2" in path.read_text().replace(" ", "")
        ]

        self.assertEqual(python2_environments, [])


if __name__ == "__main__":
    unittest.main()
