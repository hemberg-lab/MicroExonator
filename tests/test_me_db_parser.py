"""Regression tests for the mixed-format ME_DB reader in Get_annotated_microexons.py."""

import csv
import importlib.util
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
PIPELINE_MODULES = ("Bio", "pybedtools", "pyBigWig")
PIPELINE_AVAILABLE = all(
    importlib.util.find_spec(module) is not None for module in PIPELINE_MODULES
)


def _genome_sequence():
    """A 1,300 nt chr1 with AG/GT splice sites around three database microexons."""
    rng = random.Random(7)
    sequence = [rng.choice("ACT") for _ in range(1300)]  # no G: no accidental GT/AG

    def place(position, bases):
        sequence[position:position + len(bases)] = list(bases)

    place(498, "ag")   # BED6 microexon 500-512, soft-masked acceptor
    place(512, "gt")   # and soft-masked donor
    place(698, "AG")   # BED12 microexon 700-712
    place(712, "GT")
    return "".join(sequence)


@unittest.skipUnless(
    PIPELINE_AVAILABLE, "requires a pipeline environment with Biopython, pybedtools and pyBigWig"
)
class MixedFormatMicroexonDatabaseTests(unittest.TestCase):
    def run_script(self, database_rows):
        with tempfile.TemporaryDirectory() as directory:
            workdir = pathlib.Path(directory)
            (workdir / "data").mkdir()
            genome = workdir / "genome.fa"
            genome.write_text(">chr1\n" + _genome_sequence() + "\n")
            annotation = workdir / "annotation.bed12"
            # One annotated transcript: exons 100-200 and 1000-1100, intron 200-1000.
            annotation.write_text(
                "chr1\t100\t1100\ttx1\t0\t+\t100\t1100\t0\t2\t100,100\t0,900\n"
            )
            database = workdir / "ME_DB.txt"
            database.write_text(
                "".join("\t".join(map(str, row)) + "\n" for row in database_rows)
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY / "src" / "Get_annotated_microexons.py"),
                    str(genome),
                    "NA",
                    str(annotation),
                    str(REPOSITORY / "PWM" / "Mouse" / "mm10_GT_AG_U2_5.good.matrix"),
                    str(REPOSITORY / "PWM" / "Mouse" / "mm10_GT_AG_U2_3.good.matrix"),
                    "NA",
                    "30",
                    str(database),
                ],
                cwd=workdir,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with open(workdir / "data" / "DB.ME_centric") as handle:
                return {row[0] for row in csv.reader(handle, delimiter="\t")}

    def test_bed6_row_with_soft_masked_splice_sites_is_accepted(self):
        microexons = self.run_script([("chr1", 500, 512, "me_bed6", 0, "+")])

        self.assertEqual(microexons, {"chr1_+_500_512"})

    def test_one_column_row_does_not_change_length_cap_for_later_rows(self):
        """A 3 nt ID row must not turn ME_len into 3 for the BED12 row after it."""
        microexons = self.run_script(
            [
                ("chr1_+_600_603",),
                ("chr1", 100, 1100, "db_tx", 0, "+", 100, 1100, 0, 3, "100,12,100", "0,600,900"),
            ]
        )

        self.assertEqual(microexons, {"chr1_+_600_603", "chr1_+_700_712"})

    def test_all_three_formats_in_one_database(self):
        microexons = self.run_script(
            [
                ("chr1", 500, 512, "me_bed6", 0, "+"),
                ("chr1_+_600_603",),
                ("chr1", 100, 1100, "db_tx", 0, "+", 100, 1100, 0, 3, "100,12,100", "0,600,900"),
            ]
        )

        self.assertEqual(
            microexons, {"chr1_+_500_512", "chr1_+_600_603", "chr1_+_700_712"}
        )


if __name__ == "__main__":
    unittest.main()
