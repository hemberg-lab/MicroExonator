import tempfile
import unittest
from pathlib import Path

from src.insert_microexons_gtf import insert, read_me_db


def gtf_line(feature, start, end, transcript, strand="+", tags=()):
    attributes = 'gene_id "g1"; transcript_id "{}";'.format(transcript)
    attributes += "".join(' tag "{}";'.format(tag) for tag in tags)
    return "chr1\tsrc\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(feature, start, end, strand, attributes)


class ReadMeDbTests(unittest.TestCase):
    def test_bed12_bed6_and_id_lines_give_1_based_microexons(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "me_db.txt"
            path.write_text(
                "chr1\t99\t400\tME\t0\t+\t99\t400\t0\t3\t101,12,100,\t0,150,201,\n"
                "chr2\t49\t55\tx\t0\t-\n"
                "chr3_+_9_18\n")
            self.assertEqual(read_me_db(path), {("chr1", "+", 250, 261),
                                                ("chr2", "-", 50, 55), ("chr3", "+", 10, 18)})


class InsertTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.gtf = Path(self.directory.name) / "base.gtf"
        self.gtf.write_text(
            gtf_line("transcript", 100, 600, "t_long", tags=("basic",)) +
            gtf_line("exon", 100, 200, "t_long", tags=("basic",)) +
            gtf_line("exon", 301, 400, "t_long", tags=("basic",)) +
            gtf_line("exon", 501, 600, "t_long", tags=("basic",)) +
            gtf_line("transcript", 100, 400, "t_mane", tags=("MANE_Select",)) +
            gtf_line("exon", 100, 200, "t_mane", tags=("MANE_Select",)) +
            gtf_line("CDS", 150, 200, "t_mane", tags=("MANE_Select",)) +
            gtf_line("exon", 301, 400, "t_mane", tags=("MANE_Select",)))

    def tearDown(self):
        self.directory.cleanup()

    def test_microexon_goes_into_one_ranked_host_and_the_base_is_kept(self):
        base = self.gtf.read_text()
        text, report = insert(self.gtf, {("chr1", "+", 250, 261)})
        self.assertTrue(text.startswith(base))
        added = text[len(base):].splitlines()
        # MANE wins over the basic transcript; CDS lines are not copied
        self.assertTrue(all('transcript_id "t_mane_ME_250_261"' in line for line in added))
        self.assertEqual([line.split("\t")[2] for line in added], ["transcript", "exon", "exon", "exon"])
        self.assertEqual([(int(line.split("\t")[3]), int(line.split("\t")[4])) for line in added[1:]],
                         [(100, 200), (250, 261), (301, 400)])
        self.assertIn('tag "microexon_insert";', added[1])
        self.assertIn("inserted\tt_mane", report)

    def test_present_and_hostless_microexons_are_reported_not_added(self):
        base = self.gtf.read_text()
        text, report = insert(self.gtf, {("chr1", "+", 301, 400), ("chr1", "+", 700, 710),
                                         ("chr1", "-", 250, 261)})
        self.assertEqual(text, base)
        self.assertIn("301\t400\tpresent", report)
        self.assertIn("700\t710\tno_host", report)
        self.assertIn("-\t250\t261\tno_host", report)

    def test_microexon_touching_an_intron_end_has_no_host(self):
        _, report = insert(self.gtf, {("chr1", "+", 201, 210)})
        self.assertIn("no_host", report)


if __name__ == "__main__":
    unittest.main()
