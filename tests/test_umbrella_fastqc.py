"""Kept per-run QC summaries, without running FastQC itself."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src.umbrella_fastqc import fastqc, keep


class FastqcTests(unittest.TestCase):
    def test_subsamples_each_mate_and_keeps_only_fastqc_data(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reads = []
            for mate in (1, 2):
                path = root / "R{}.fastq.gz".format(mate)
                with gzip.open(path, "wt") as stream:
                    for number in range(5):
                        stream.write("@r{}\nACGT\n+\nIIII\n".format(number))
                reads.append(str(path))
            seen = {}

            def fake_fastqc(command, check):
                work = Path(command[command.index("-o") + 1])
                for fastq in command[command.index("-o") + 2:]:
                    name = Path(fastq).stem
                    seen[name] = Path(fastq).read_text().count("\n") // 4
                    (work / (name + "_fastqc")).mkdir()
                    (work / (name + "_fastqc") / "fastqc_data.txt").write_text("##FastQC\t0.12.1\n")

            with patch("src.umbrella_fastqc.subprocess.run", side_effect=fake_fastqc):
                fastqc(reads, "SRR1", root / "out", reads_per_file=2, scratch=root)
            self.assertEqual(seen, {"SRR1_R1": 2, "SRR1_R2": 2})
            self.assertEqual(sorted(p.relative_to(root / "out").as_posix()
                                    for p in (root / "out").rglob("*") if p.is_file()),
                             ["SRR1_R1/fastqc_data.txt", "SRR1_R2/fastqc_data.txt"])
            self.assertFalse(list(root.glob("fastqc-*")))


class KeepTests(unittest.TestCase):
    def test_copies_summaries_in_multiqc_layout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "hisat2.txt").write_text("HISAT2 summary stats:\n")
            (root / "fc.summary").write_text("Status\tumbrella/work/x/align/aligned.bam\nAssigned\t10\n")
            (root / "salmon" / "aux_info").mkdir(parents=True)
            (root / "salmon" / "aux_info" / "meta_info.json").write_text(json.dumps({"percent_mapped": 90}))
            (root / "salmon" / "lib_format_counts.json").write_text("{}")
            (root / "quant.map.gz").write_bytes(gzip.compress(b"map\n"))
            out = root / "kept"
            keep("SRR1", root / "hisat2.txt", root / "fc.summary", root / "salmon",
                 root / "quant.map.gz", out)
            self.assertEqual((out / "SRR1.featurecounts.txt.summary").read_text(),
                             "Status\tSRR1\nAssigned\t10\n")
            for name in ("SRR1.hisat2.summary.txt", "salmon/SRR1/aux_info/meta_info.json",
                         "salmon/SRR1/lib_format_counts.json", "SRR1.whippet.map.gz"):
                self.assertTrue((out / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
