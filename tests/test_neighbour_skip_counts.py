"""Tests for exclusion counts of consecutive microexons (neighbour skips)."""

import csv
import importlib.util
import pathlib
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import neighbour_skip_counts as nsc  # noqa: E402

BIOPYTHON_AVAILABLE = importlib.util.find_spec("Bio") is not None

# E1 ends at 1000, m1 = 1200-1209, m2 = 1400-1418, E2 starts at 2000 (+ strand).
M1 = "chr1_+_1200_1209"
M2 = "chr1_+_1400_1418"
M1_SEQ = "AATCTATTG"
M2_SEQ = "GTGGTACACCAATAAGAG"


class NeighbourSkipHelperTests(unittest.TestCase):
    def test_parsers_accept_underscores_in_chromosome_names(self):
        self.assertEqual(nsc.parse_intron("chrUn_GL1:5-90"), ("chrUn_GL1", "-", 5, 90))
        self.assertEqual(nsc.parse_microexon("chrUn_GL1_-_10_19"), ("chrUn_GL1", "-", 10, 19))

    def test_index_returns_only_microexons_fully_inside_the_interval(self):
        index = nsc.MicroexonIndex([M1, M2, "chr1_+_990_1003", "chr1_-_1200_1209"])
        self.assertEqual(index.inside("chr1", "+", 1000, 1400), [M1])
        self.assertEqual(index.inside("chr1", "+", 1000, 2000), [M1, M2])
        self.assertEqual(index.inside("chr2", "+", 0, 5000), [])

    def test_plus_strand_upstream_junction_skips_the_previous_microexon(self):
        index = nsc.MicroexonIndex([M1, M2])
        # Read on m2's E1|m2|E2 tag crossing E1|m2 with 20 nt anchors.
        skipped = nsc.skipped_by_read(index, "chr1:1000+2000", M2, 80, 40, 100, 18)
        self.assertEqual(skipped, {M1: "chr1:1000+1400"})

    def test_read_inside_the_microexon_side_only_crosses_nothing(self):
        index = nsc.MicroexonIndex([M1, M2])
        self.assertEqual(nsc.skipped_by_read(index, "chr1:1000+2000", M2, 95, 10, 100, 18), {})

    def test_junction_between_adjacent_microexons_skips_nothing(self):
        index = nsc.MicroexonIndex([M1, M2])
        # m2's m1|m2|E2 tag: its upstream junction is m1|m2.
        self.assertEqual(nsc.skipped_by_read(index, "chr1:1209+2000", M2, 80, 40, 100, 18), {})

    def test_minus_strand_uses_transcript_order(self):
        # Minus strand: transcript order is E1 (>2000), m2, m1, E2 (<1000).
        index = nsc.MicroexonIndex(["chr1_-_1200_1209", "chr1_-_1400_1418"])
        # Read on m1's E1|m1|E2 tag crossing its upstream junction E1|m1 (1209-2000).
        up = nsc.skipped_by_read(index, "chr1:1000-2000", "chr1_-_1200_1209", 80, 40, 100, 9)
        self.assertEqual(up, {"chr1_-_1400_1418": "chr1:1209-2000"})
        # Read on m2's tag crossing its downstream junction m2|E2 (1000-1400).
        down = nsc.skipped_by_read(index, "chr1:1000-2000", "chr1_-_1400_1418", 110, 40, 100, 18)
        self.assertEqual(down, {"chr1_-_1200_1209": "chr1:1000-1400"})

    def test_run_of_three_counts_each_skipped_microexon_once(self):
        m3 = "chr1_+_1600_1612"
        index = nsc.MicroexonIndex([M1, M2, m3])
        # Read on m3's E1|m3|E2 tag crossing E1|m3 and m3|E2 skips m1 and m2 once each.
        skipped = nsc.skipped_by_read(index, "chr1:1000+2000", m3, 85, 40, 100, 12)
        self.assertEqual(skipped, {M1: "chr1:1000+1600", M2: "chr1:1000+1600"})


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "needs the pipeline Biopython environment")
class MicroexonCoverageScriptTests(unittest.TestCase):
    def run_script(self, m1_introns, neighbour_skips=None):
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            tags = directory / "tags.fa"
            tags.write_text(
                ">chr1:1000+1400|tx1|100_{}_100\nACGT\n".format(M1_SEQ)
                + ">chr1:1000+2000|tx1|100_{}_100\nACGT\n".format(M1_SEQ)
                + ">chr1:1209+2000|tx1|100_{}_100\nACGT\n".format(M2_SEQ)
                + ">chr1:1000+2000|tx1|100_{}_100\nACGT\n".format(M2_SEQ)
                + ">chr1:1000+2000|tx1|100_100\nACGT\n"
            )
            centric = directory / "TOTAL.ME_centric.txt"
            rows = [
                (M1, m1_introns, M1_SEQ, M1 + "|80|0"),
                (M2, "chr1:1209+2000,chr1:1000+2000", M2_SEQ, M2 + "|80|0"),
            ]
            centric.write_text(
                "".join(
                    "\t".join([me, "tx1", "0", sjs, "0", str(len(seq)), seq, "1", "80", "0", "0", total])
                    + "\n"
                    for me, sjs, seq, total in rows
                )
            )
            bed12 = directory / "genes.bed12"
            bed12.write_text("chr1\t500\t2500\ttx1\t0\t+\t500\t2500\t0\t2\t500,500,\t0,1500,\n")
            reads = directory / "S1.sam.pre_processed.filter1"
            reads.write_text(
                "".join(
                    "\t".join([name, "0", tag, str(start), "40M", "A" * 40, "I" * 40]) + "\n"
                    for name, tag, start in [
                        # E1|m2 read on m2's E1|m2|E2 tag: m2 inclusion, m1 exclusion.
                        ("r1", "chr1:1000+2000|tx1|100_{}_100".format(M2_SEQ), 80),
                        # m1|E2 read on m1's E1|m1|E2 tag: m1 inclusion, m2 exclusion.
                        ("r2", "chr1:1000+2000|tx1|100_{}_100".format(M1_SEQ), 90),
                        # m1|m2 read on m1's E1|m1|m2 tag: skips nothing.
                        ("r3", "chr1:1000+1400|tx1|100_{}_100".format(M1_SEQ), 90),
                        # Double skip E1|E2: exclusion for both.
                        ("r4", "chr1:1000+2000|tx1|100_100", 80),
                    ]
                )
            )
            command = [
                sys.executable,
                str(REPOSITORY / "src" / "ME_SJ_coverage.py"),
                str(tags), str(centric), str(bed12), str(reads), "30",
            ]
            if neighbour_skips is not None:
                command.append(neighbour_skips)
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            output = {}
            for row in csv.reader(result.stdout.splitlines(), delimiter="\t"):
                output[row[1]] = {
                    "SJs": row[2].split(","),
                    "ME_covs": row[3].split(","),
                    "sum_ME": int(row[4]),
                    "SJ_covs": list(map(int, row[8].split(","))),
                    "sum_SJ": int(row[9]),
                }
            return output

    def test_single_skip_reads_count_as_exclusion_for_the_neighbour(self):
        output = self.run_script("chr1:1000+1400,chr1:1000+2000")
        self.assertEqual(output[M1]["SJ_covs"], [1, 1])
        self.assertEqual(output[M2]["SJ_covs"], [1, 1])
        self.assertEqual(output[M1]["sum_ME"], 2)
        self.assertEqual(output[M2]["sum_ME"], 1)

    def test_option_off_reproduces_the_old_counts(self):
        output = self.run_script("chr1:1000+1400,chr1:1000+2000", "F")
        self.assertEqual(output[M1]["SJ_covs"], [0, 1])
        self.assertEqual(output[M2]["SJ_covs"], [0, 1])

    def test_skipping_junction_missing_from_the_event_is_appended(self):
        output = self.run_script("chr1:1000+2000")
        self.assertEqual(output[M1]["SJs"], ["chr1:1000+2000", "chr1:1000+1400"])
        self.assertEqual(output[M1]["ME_covs"], ["1", "0"])
        self.assertEqual(output[M1]["SJ_covs"], [1, 1])
        self.assertEqual(output[M1]["sum_SJ"], 2)


if __name__ == "__main__":
    unittest.main()
