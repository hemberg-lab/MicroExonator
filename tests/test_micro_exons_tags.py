"""Tests for building inclusion tags of discovered microexons (Micro_exons_tags.py)."""

import importlib.util
import pathlib
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
BIOPYTHON_AVAILABLE = importlib.util.find_spec("Bio") is not None


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "needs the pipeline Biopython environment")
class MicroExonTagTests(unittest.TestCase):
    def test_every_skipping_tag_path_of_the_intron_gets_an_inclusion_tag(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            tags = directory / "ME_TAGs.fa"
            tags.write_text(
                ">chr1:200+1000|txA|4_4\nAAAACCCC\n"
                ">chr1:200+1000|txB|4_4\nGGGGCCCC\n"
                ">chr2:5+900|txC|4_4\nTTTTAAAA\n"
            )
            centric = directory / "ME_centric"
            centric.write_text("\t".join(
                ["chr1_+_500_503", "txA", "0", "chr1:200+1000", "0", "3", "TAG", "1", "80", "0", "0", "chr1_+_500_503|80|0"]
            ) + "\n")
            result = subprocess.run(
                [sys.executable, str(REPOSITORY / "src" / "Micro_exons_tags.py"), str(tags), str(centric)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = result.stdout.splitlines()
            records = dict(zip(lines[0::2], lines[1::2]))
            self.assertEqual(records[">chr1:200+1000|txA|4_TAG_4"], "AAAATAGCCCC")
            self.assertEqual(records[">chr1:200+1000|txB|4_TAG_4"], "GGGGTAGCCCC")
            self.assertIn(">chr2:5+900|txC|4_4", records)


if __name__ == "__main__":
    unittest.main()
