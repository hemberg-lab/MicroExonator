import gzip
import json
import tempfile
import unittest
from pathlib import Path

from src import circrna_catalogue as catalogue
from tests import circrna_fixture as fixture


class CatalogueTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.paths = fixture.build(self.dir / "fx")

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, **options):
        return catalogue.build(self.paths["catalogue"], self.paths["genome"], self.dir / "out", **options)

    def test_coordinates_become_ciriquant_one_based_ids(self):
        report = self.build()
        bed = (self.dir / "out" / "catalogue.bed").read_text().splitlines()
        self.assertEqual(bed, ["chr1\t201\t1150\tchr1:201|1150\t.\t+",
                               "chr1\t1001\t1800\tchr1:1001|1800\t.\t+",
                               "chr2\t2201\t2400\tchr2:2201|2400\t.\t-"])
        self.assertEqual(report["kept"], 3)
        self.assertEqual(report["source_rows"], 8)

    def test_drops_are_counted_by_reason(self):
        report = self.build()
        self.assertEqual(report["dropped"], {"contig_not_in_genome": 1, "outside_contig_or_empty": 1,
                                             "same_bsj_on_both_strands": 1})

    def test_names_of_one_bsj_are_merged(self):
        self.build()
        with gzip.open(self.dir / "out" / "annotation.tsv.gz", "rt") as stream:
            rows = {line.split("\t")[0]: line.rstrip("\n").split("\t") for line in stream if not line.startswith("circ_id")}
        self.assertEqual(rows["chr1:1001|1800"][6], "circAtlas_C1,circBase_C1")

    def test_true_circrnas_have_canonical_splice_motifs_on_both_strands(self):
        self.build()
        with gzip.open(self.dir / "out" / "annotation.tsv.gz", "rt") as stream:
            motifs = {line.split("\t")[0]: line.rstrip("\n").split("\t")[7] for line in stream}
        self.assertEqual(motifs["chr1:1001|1800"], "GT-AG")
        self.assertEqual(motifs["chr2:2201|2400"], "GT-AG")

    def test_sequence_across_the_wrap(self):
        # the pseudo-circular sequence's junction joins the circle's last base to its first
        genome = catalogue.Fasta(self.paths["genome"])
        span = genome.fetch("chr1", 1000, 1800)
        self.assertEqual(span, self.paths["chroms"]["chr1"][1000:1800])
        junction = (span * 2)[len(span) - 10:len(span) + 10]
        self.assertEqual(junction, self.paths["chroms"]["chr1"][1790:1800] + self.paths["chroms"]["chr1"][1000:1010])

    def test_span_cap_and_size_refusal(self):
        report = self.build(max_span=900)
        self.assertEqual(report["dropped"]["span_above_max_span"], 1)
        with self.assertRaisesRegex(SystemExit, "pseudo-circular reference would hold"):
            self.build(max_pseudo_bases=100)
        refused = json.loads((self.dir / "out" / "report.json").read_text())
        self.assertIn("refused", refused)

    def test_identity_is_deterministic_and_content_bound(self):
        first = self.build()["catalogue_id"]
        self.assertEqual(self.build()["catalogue_id"], first)
        self.assertNotEqual(self.build(max_span=900)["catalogue_id"], first)

    def test_read_catalogue(self):
        self.build()
        entries = catalogue.read_catalogue(self.dir / "out" / "catalogue.bed")
        self.assertEqual(entries["chr2:2201|2400"], ("chr2", 2201, 2400, "-"))


CHAIN = """chain 1000 chrA 5000 + 0 1000 chrA 6000 + 100 1090 1
400 10 0
590

chain 900 chrB 3000 + 0 500 chrB 2000 - 0 500 2
500

chain 800 chrC 1000 + 100 300 chrC 1000 + 100 300 3
200

chain 700 chrC 1000 + 100 200 chrX 1000 + 500 600 4
100
"""


class LiftTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        with gzip.open(self.dir / "test.over.chain.gz", "wt") as stream:
            stream.write(CHAIN)
        (self.dir / "hg19.bed").write_text("".join([
            "chrA\t50\t300\tsame_span\t0\t+\n",
            "chrA\t300\t600\tspan_shrinks\t0\t+\n",
            "chrA\t395\t405\tend_in_gap\t0\t+\n",
            "chrB\t100\t200\treverse_chain\t0\t+\n",
            "chrC\t150\t160\tambiguous\t0\t-\n",
            "chrZ\t1\t50\tno_chain\t0\t+\n"]))
        self.report = catalogue.lift(self.dir / "hg19.bed", self.dir / "test.over.chain.gz",
                                     self.dir / "hg38.bed", self.dir / "failed.tsv")

    def tearDown(self):
        self.tmp.cleanup()

    def test_boundaries_lift_one_by_one(self):
        rows = [line.split("\t") for line in (self.dir / "hg38.bed").read_text().splitlines()]
        self.assertEqual(rows[0], ["chrA", "150", "400", "same_span", "0", "+"])
        # 599 sits after a 10-base deletion in the target: the circle loses 10 bases
        self.assertEqual(rows[1], ["chrA", "400", "690", "span_shrinks", "0", "+"])
        # a reverse-strand chain flips both the order of the boundaries and the strand
        self.assertEqual(rows[2], ["chrB", "1800", "1900", "reverse_chain", "0", "-"])
        self.assertEqual(self.report["lifted"], 3)
        self.assertEqual(self.report["span_changed"], 1)

    def test_failures_are_kept_with_reasons(self):
        self.assertEqual(self.report["failed"], {"boundary_in_two_chains": 1, "boundary_not_in_aligned_block": 1,
                                                 "contig_not_in_chain": 1})
        failed = (self.dir / "failed.tsv").read_text().splitlines()
        self.assertEqual(len(failed), 4)
        self.assertIn("end_in_gap", failed[1])

    def test_lifted_rows_feed_the_build(self):
        # a lifted source and a native hg38 source build into one catalogue
        fx = fixture.build(self.dir / "fx")
        (self.dir / "lifted.bed").write_text("chr1\t1000\t1800\tcircBase_C1_lifted\t0\t+\n")
        report = catalogue.build([fx["catalogue"], self.dir / "lifted.bed"], fx["genome"], self.dir / "out")
        self.assertEqual(report["source_rows"], 9)
        with gzip.open(self.dir / "out" / "annotation.tsv.gz", "rt") as stream:
            names = {line.split("\t")[0]: line.split("\t")[6] for line in stream}
        self.assertEqual(names["chr1:1001|1800"], "circAtlas_C1,circBase_C1,circBase_C1_lifted")


class ConvertTests(unittest.TestCase):
    """Each database's own layout and coordinate convention (checked against hg38/hg19
    splice motifs on 25 sampled rows each, 25 October 2026)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def convert(self, text, fmt, tissue=None):
        source = self.dir / ("in." + fmt)
        source.write_text(text)
        report = catalogue.convert(source, self.dir / "out.bed", fmt, tissue)
        return [line.split("\t") for line in (self.dir / "out.bed").read_text().splitlines()], report

    def test_circatlas_start_is_one_based(self):
        rows, _ = self.convert("Chro\tStart\tEnd\tStrand\tcircAltas_ID\n"
                               "chr1\t926703\t927095\t-\thsa-SAMD11_0013\n", "circatlas")
        self.assertEqual(rows, [["chr1", "926702", "927095", "hsa-SAMD11_0013", "0", "-"]])

    def test_circbase_start_is_zero_based(self):
        rows, _ = self.convert("# chrom\tstart\tend\tstrand\tcircRNA ID\tgenomic length\n"
                               "chr10\t28196591\t28233891\t-\thsa_circ_0018046\t37300\n", "circbase")
        self.assertEqual(rows, [["chr10", "28196591", "28233891", "hsa_circ_0018046", "0", "-"]])

    def test_circpedia_location_and_tissue_filter(self):
        text = ("circID\tLocation\tcircname\tTGS\tspecies\n"
                "CIRCHSA_A_1\tchr22:20400631-20401372(+)\tcircA\tCentral nervous system: Brain[isoCirc]\tHomo-sapiens\n"
                "CIRCHSA_B_1\tchr2:97584943-97585398(-)\tcircB\tnone\tHomo-sapiens\n")
        rows, report = self.convert(text, "circpedia", "Central nervous system")
        self.assertEqual(rows, [["chr22", "20400631", "20401372", "CIRCHSA_A_1", "0", "+"]])
        self.assertEqual((report["kept"], report["skipped"]), (1, 1))
        rows, _ = self.convert(text, "circpedia")
        self.assertEqual(len(rows), 2)

    def test_annotation_sources_only_name_existing_bsjs(self):
        fx = fixture.build(self.dir / "fx")
        (self.dir / "names.bed").write_text("chr1\t1000\t1800\tcircBase_extra\t0\t+\n"
                                            "chr1\t1000\t1800\twrong_strand\t0\t-\n"
                                            "chr1\t2000\t2500\tnot_in_core\t0\t+\n")
        report = catalogue.build(fx["catalogue"], fx["genome"], self.dir / "out", annotate=[self.dir / "names.bed"])
        self.assertEqual((report["kept"], report["annotation_names"]), (3, 1))
        with gzip.open(self.dir / "out" / "annotation.tsv.gz", "rt") as stream:
            names = {line.split("\t")[0]: line.split("\t")[6] for line in stream}
        self.assertIn("circBase_extra", names["chr1:1001|1800"])
        self.assertNotIn("wrong_strand", names["chr1:1001|1800"])


if __name__ == "__main__":
    unittest.main()
