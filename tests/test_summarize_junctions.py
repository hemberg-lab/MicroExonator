import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src.summarize_junctions import (Catalog, JunctionCounter, build_catalog,
                                     capture_text, infer_strandedness,
                                     read_junctions, run_capture_rates,
                                     shard_text, strandedness_argument)


def sam(qname, flag, pos, cigar, mate_pos=0, nh=1, xs="+", chrom="chr1", mate_chrom="="):
    tags = ["NH:i:{}".format(nh)] + (["XS:A:{}".format(xs)] if xs else [])
    return "\t".join([qname, str(flag), chrom, str(pos), "60", cigar, mate_chrom,
                      str(mate_pos), "0", "A", "I"] + tags) + "\n"


class ReadJunctionTests(unittest.TestCase):
    def test_anchors_stop_at_neighbouring_introns(self):
        # 10M 100N 5M 50N 20M from 1000: introns 1010-1109 and 1115-1164.
        self.assertEqual(read_junctions(1000, "10M100N5M50N20M"),
                         [(1010, 1109, 5), (1115, 1164, 5)])

    def test_deletions_advance_reference_insertions_and_clips_do_not(self):
        self.assertEqual(read_junctions(1000, "2S10M2D3I10M100N12M"),
                         [(1022, 1121, 12)])


class JunctionCounterTests(unittest.TestCase):
    def count(self, lines, min_anchor=8):
        counter = JunctionCounter(min_anchor)
        for line in lines:
            counter.add(line)
        return counter.finish()

    def test_one_count_per_fragment_and_junction(self):
        counter = self.count([
            sam("f1", 0x1 | 0x40, 1000, "20M100N30M", mate_pos=1005),
            sam("f1", 0x1 | 0x80 | 0x10, 1005, "15M100N35M", mate_pos=1000),
        ])
        self.assertEqual(counter.table[("chr1", 1020, 1119, "+")], [1, 0, 0, 20])

    def test_unique_multimapping_and_short_anchor_are_separate(self):
        counter = self.count([
            sam("u", 0, 1000, "20M100N30M"),
            sam("m", 0, 1000, "20M100N30M", nh=3),
            sam("s", 0, 1013, "7M100N40M"),
        ])
        self.assertEqual(counter.table[("chr1", 1020, 1119, "+")], [1, 1, 1, 20])

    def test_secondary_supplementary_qcfail_and_unmapped_are_skipped(self):
        counter = self.count([sam("x", flag, 1000, "20M100N30M")
                              for flag in (0x4, 0x100, 0x200, 0x800)])
        self.assertEqual(dict(counter.table), {})

    def test_reads_without_xs_are_ambiguous_and_excluded(self):
        counter = self.count([sam("a", 0, 1000, "20M100N30M", xs=None)])
        self.assertEqual(dict(counter.table), {})
        self.assertEqual(counter.summary["ambiguous_strand"], 1)

    def test_mate_filtered_out_still_counts_once(self):
        counter = self.count([sam("f", 0x1 | 0x40, 1000, "20M100N30M", mate_pos=1300)])
        self.assertEqual(counter.table[("chr1", 1020, 1119, "+")], [1, 0, 0, 20])

    def test_strandedness_from_xs_with_read2_flipped(self):
        # dUTP (firststrand): read 1 antisense to the transcript.
        lines = []
        for i in range(600):
            lines.append(sam("r{}".format(i), 0x1 | 0x40 | 0x10, 1000, "20M100N30M",
                             mate_pos=1400, mate_chrom="chr2"))
            lines.append(sam("q{}".format(i), 0x1 | 0x80, 2000, "20M100N30M",
                             mate_pos=1400, mate_chrom="chr2"))
        counter = self.count(lines)
        self.assertEqual(counter.summary["sense"], 0)
        self.assertEqual(counter.summary["inferred_strandedness"], "firststrand")

    def test_strandedness_needs_enough_evidence(self):
        self.assertEqual(infer_strandedness(10, 0), "undetermined")
        self.assertEqual(infer_strandedness(500, 500), "unstranded")
        self.assertEqual(infer_strandedness(900, 100), "secondstrand")


class StrandArgumentTests(unittest.TestCase):
    def test_declared_wins_and_auto_uses_inferred(self):
        self.assertEqual(strandedness_argument("featurecounts", "firststrand", None), "2")
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / "s.json"
            summary.write_text(json.dumps({"inferred_strandedness": "secondstrand"}))
            self.assertEqual(strandedness_argument("featurecounts", "auto", summary), "1")
            summary.write_text(json.dumps({"inferred_strandedness": "undetermined"}))
            self.assertEqual(strandedness_argument("featurecounts", "auto", summary), "0")
        self.assertEqual(strandedness_argument("hisat2", "firststrand", None, "PE"), "RR")


class CatalogShardCaptureTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        root = Path(self.directory.name)
        # Whippet index: exons 100-200, 301-400 (+) -> intron 201-300.
        #                exons 100-200, 501-600      -> intron 201-500.
        (root / "whippet.gtf").write_text(
            'chr1\tx\texon\t100\t200\t.\t+\t.\ttranscript_id "a";\n'
            'chr1\tx\texon\t301\t400\t.\t+\t.\ttranscript_id "a";\n'
            'chr1\tx\texon\t100\t180\t.\t+\t.\ttranscript_id "b";\n'
            'chr1\tx\texon\t501\t600\t.\t+\t.\ttranscript_id "b";\n')
        # Full annotation adds a pruned acceptor: intron 201-303.
        (root / "annotation.gtf").write_text(
            'chr1\tx\texon\t100\t200\t.\t+\t.\ttranscript_id "c";\n'
            'chr1\tx\texon\t304\t400\t.\t+\t.\ttranscript_id "c";\n')
        # Hint only: intron 201-700 (0-based last exon base 199, first exon base 700).
        (root / "sites.txt").write_text("chr1\t199\t700\t+\n")
        catalog = root / "catalog.tsv.gz"
        catalog.write_bytes(gzip.compress(build_catalog(
            root / "whippet.gtf", root / "annotation.gtf", root / "sites.txt").encode()))
        self.catalog = Catalog(catalog)
        self.root = root

    def tearDown(self):
        self.directory.cleanup()

    def test_labels_include_novel_combinations_of_index_sites(self):
        self.assertEqual(self.catalog.label(("chr1", 201, 300, "+")), "index_intron")
        self.assertEqual(self.catalog.label(("chr1", 181, 300, "+")), "index_sites")
        self.assertEqual(self.catalog.label(("chr1", 201, 303, "+")), "annotation")
        self.assertEqual(self.catalog.label(("chr1", 201, 700, "+")), "hint")
        self.assertEqual(self.catalog.label(("chr1", 201, 800, "+")), "novel")
        self.assertEqual(self.catalog.label(("chr1", 201, 300, "-")), "novel")

    def write_run(self, name, rows):
        path = self.root / (name + ".tsv")
        path.write_text("chrom\tstart\tend\tstrand\tunique\tmulti\tshort_anchor\tmax_anchor\n" +
                        "".join("\t".join(map(str, row)) + "\n" for row in rows))
        return str(path)

    def test_shard_and_capture_flag_the_pruned_acceptor(self):
        tables = {
            "run2": self.write_run("run2", [("chr1", 201, 300, "+", 2, 0, 0, 30)]),
            "run1": self.write_run("run1", [("chr1", 201, 300, "+", 6, 1, 0, 40),
                                            ("chr1", 201, 303, "+", 4, 0, 2, 25)]),
        }
        shard = self.root / "shard.tsv"
        shard.write_text(shard_text(tables, self.catalog))
        self.assertEqual(shard.read_text().splitlines(), [
            "chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\tshort_anchor_total\trun1\trun2",
            "chr1\t201\t300\t+\tindex_intron\t40\t1\t0\t6\t2",
            "chr1\t201\t303\t+\tannotation\t25\t0\t2\t4\t0",
        ])
        capture = capture_text(shard, self.catalog).splitlines()
        self.assertEqual(len(capture), 2)
        self.assertEqual(capture[1].split("\t"), [
            "chr1", "201", "+", "left", "chr1:201-303:+", "4",
            "6", "4", "0.6", "2", "0", "1"])
        self.assertEqual(run_capture_rates(shard), {"run1": 0.6, "run2": 1.0})

    def test_cli_entrypoint_loads_without_pythonpath(self):
        script = Path(__file__).resolve().parents[1] / "src/summarize_junctions.py"
        result = subprocess.run([sys.executable, str(script), "--help"], text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(result.returncode, 0, result.stdout)


if __name__ == "__main__":
    unittest.main()
