"""Tests for synthetic skipping tags: pure helpers and the annotation script."""

import csv
import importlib.util
import pathlib
import random
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

import synthetic_skip_tags as sst  # noqa: E402

PIPELINE_AVAILABLE = all(
    importlib.util.find_spec(module) is not None
    for module in ("Bio", "pybedtools", "pyBigWig")
)


class SyntheticSkipTagHelperTests(unittest.TestCase):
    def test_existing_tag_introns_are_read_from_two_anchor_headers_only(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = pathlib.Path(directory) / "tags.fa"
            fasta.write_text(
                ">chr1:200+1000|tx1|100_100\nACGT\n"
                ">chr1:200+1000|tx1|100_ACG_100\nACGT\n"
                ">chr2:5-9|tx2|40_60\nACGT\n"
            )
            self.assertEqual(
                sst.load_tag_introns(fasta), {"chr1:200+1000", "chr2:5-9"}
            )

    def test_intron_without_tag_and_not_bordering_a_microexon_needs_a_tag(self):
        boundaries = sst.microexon_boundaries([("chr1", "+", 500, 512)])
        self.assertTrue(
            sst.needs_synthetic_tag("chr1", "+", 200, 1000, set(), set(), boundaries)
        )

    def test_existing_or_already_written_introns_get_no_tag(self):
        boundaries = sst.microexon_boundaries([("chr1", "+", 500, 512)])
        self.assertFalse(
            sst.needs_synthetic_tag("chr1", "+", 200, 1000, {"chr1:200+1000"}, set(), boundaries)
        )
        self.assertFalse(
            sst.needs_synthetic_tag("chr1", "+", 200, 1000, set(), {"chr1:200+1000"}, boundaries)
        )

    def test_introns_bordering_another_microexon_get_no_tag(self):
        """Consecutive microexons: the skip junction is the neighbour's inclusion junction."""
        boundaries = sst.microexon_boundaries(
            [("chr1", "+", 500, 512), ("chr1", "+", 600, 609)]
        )
        # m1's containing intron ends where m2 starts
        self.assertFalse(sst.needs_synthetic_tag("chr1", "+", 200, 600, set(), set(), boundaries))
        # m2's containing intron starts where m1 ends
        self.assertFalse(sst.needs_synthetic_tag("chr1", "+", 512, 1000, set(), set(), boundaries))
        # the other strand is unaffected
        self.assertTrue(sst.needs_synthetic_tag("chr1", "-", 200, 600, set(), set(), boundaries))

    def test_plus_strand_record_joins_flanks_in_genomic_order(self):
        header, sequence = sst.skip_tag_record("chr1", "+", 200, 1000, "tx1", "aacc", "GGT")
        self.assertEqual(header, "chr1:200+1000|tx1|4_3")
        self.assertEqual(sequence, "AACCGGT")

    def test_minus_strand_record_is_reverse_complemented_with_swapped_anchors(self):
        header, sequence = sst.skip_tag_record("chr1", "-", 200, 1000, "tx1", "AACC", "GGT")
        self.assertEqual(header, "chr1:200-1000|tx1|3_4")
        self.assertEqual(sequence, "ACCGGTT")

    def test_microexon_ids_with_underscored_chromosomes_are_parsed(self):
        self.assertEqual(
            sst.parse_microexon_id("chr1_GL456210_random_-_10_19"),
            ("chr1_GL456210_random", "-", 10, 19),
        )


def _genome():
    """chr1 without G except the AG/GT splice sites placed around the microexons."""
    rng = random.Random(11)
    sequence = [rng.choice("ACT") for _ in range(1300)]
    for position, bases in ((498, "AG"), (512, "GT"), (598, "AG"), (609, "GT")):
        sequence[position:position + 2] = list(bases)
    return "".join(sequence)


@unittest.skipUnless(
    PIPELINE_AVAILABLE, "requires a pipeline environment with Biopython, pybedtools and pyBigWig"
)
class SyntheticSkipTagGenerationTests(unittest.TestCase):
    """Runs Get_annotated_microexons.py on a 1.3 kb synthetic chromosome.

    Exons used: E1 100-200, m1 500-512 (12 nt), m2 600-609 (9 nt), E2 1000-1100.
    """

    def run_script(self, transcripts, strand="+", existing_tags="", enabled=True):
        with tempfile.TemporaryDirectory() as directory:
            workdir = pathlib.Path(directory)
            (workdir / "data").mkdir()
            genome = _genome()
            if strand == "-":
                # Mirror the splice-site signals so the minus-strand transcript is canonical.
                complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
                genome = list(genome)
                for position, bases in ((498, "AC"), (512, "CT"), (598, "AC"), (609, "CT")):
                    genome[position:position + 2] = list(bases)
                genome = "".join(genome)
            (workdir / "genome.fa").write_text(">chr1\n" + genome + "\n")
            with open(workdir / "annotation.bed12", "w") as handle:
                for name, exons in transcripts.items():
                    start, end = exons[0][0], exons[-1][1]
                    handle.write("\t".join(map(str, [
                        "chr1", start, end, name, 0, strand, start, end, 0, len(exons),
                        ",".join(str(e - s) for s, e in exons),
                        ",".join(str(s - start) for s, e in exons),
                    ])) + "\n")
            (workdir / "empty.bed").write_text("")
            (workdir / "ME_TAGs.fa").write_text(existing_tags)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY / "src" / "Get_annotated_microexons.py"),
                    str(workdir / "genome.fa"), "NA", str(workdir / "annotation.bed12"),
                    str(REPOSITORY / "PWM" / "Mouse" / "mm10_GT_AG_U2_5.good.matrix"),
                    str(REPOSITORY / "PWM" / "Mouse" / "mm10_GT_AG_U2_3.good.matrix"),
                    "NA", "30", str(workdir / "empty.bed"),
                    str(workdir / "ME_TAGs.fa") if enabled else "NA",
                ],
                cwd=workdir, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            headers = [
                line[1:].strip()
                for line in (workdir / "data" / "ME_canonical_SJ_tags.DB.fa").read_text().splitlines()
                if line.startswith(">")
            ]
            sequences = [
                line.strip()
                for line in (workdir / "data" / "ME_canonical_SJ_tags.DB.fa").read_text().splitlines()
                if not line.startswith(">")
            ]
            skip_tags = {
                h: s for h, s in zip(headers, sequences) if len(h.split("|")[2].split("_")) == 2
            }
            with open(workdir / "data" / "DB.ME_centric") as handle:
                centric = {row[0]: row[3] for row in csv.reader(handle, delimiter="\t")}
            return skip_tags, centric, genome

    def test_microexon_only_annotated_as_included_gets_one_skipping_tag(self):
        skip_tags, centric, genome = self.run_script(
            {"tx_inc": [(100, 200), (500, 512), (1000, 1100)]}
        )
        self.assertEqual(centric, {"chr1_+_500_512": "chr1:200+1000"})
        self.assertEqual(list(skip_tags), ["chr1:200+1000|tx_inc|100_100"])
        self.assertEqual(skip_tags["chr1:200+1000|tx_inc|100_100"], genome[100:200] + genome[1000:1100])

    def test_existing_skipping_tag_is_not_duplicated(self):
        skip_tags, centric, genome = self.run_script(
            {"tx_inc": [(100, 200), (500, 512), (1000, 1100)]},
            existing_tags=">chr1:200+1000|tx_skip|100_100\nACGT\n",
        )
        self.assertEqual(skip_tags, {})

    def test_option_off_writes_no_synthetic_tags(self):
        skip_tags, centric, genome = self.run_script(
            {"tx_inc": [(100, 200), (500, 512), (1000, 1100)]}, enabled=False
        )
        self.assertEqual(skip_tags, {})

    def test_consecutive_microexons_get_no_synthetic_tag_on_shared_introns(self):
        skip_tags, centric, genome = self.run_script(
            {"tx_both": [(100, 200), (500, 512), (600, 609), (1000, 1100)]}
        )
        self.assertEqual(set(centric), {"chr1_+_500_512", "chr1_+_600_609"})
        self.assertEqual(skip_tags, {})

    def test_consecutive_pair_still_gets_double_skip_tag_when_that_intron_is_annotated(self):
        skip_tags, centric, genome = self.run_script(
            {
                "tx_both": [(100, 200), (500, 512), (600, 609), (1000, 1100)],
                "tx_none": [(100, 200), (1000, 1100)],
            }
        )
        # E1-E2 is annotated but, with no existing tags supplied, still untagged.
        self.assertEqual([h.split("|")[0] for h in skip_tags], ["chr1:200+1000"])

    def test_minus_strand_tag_is_reverse_complemented(self):
        skip_tags, centric, genome = self.run_script(
            {"tx_inc": [(100, 200), (500, 512), (1000, 1100)]}, strand="-"
        )
        complement = str.maketrans("ACGT", "TGCA")
        expected = (genome[100:200] + genome[1000:1100]).translate(complement)[::-1]
        self.assertEqual(skip_tags, {"chr1:200-1000|tx_inc|100_100": expected})


if __name__ == "__main__":
    unittest.main()
